#!/usr/bin/env python3
"""
pipeline_status.py
------------------
Suivi d'un pipeline à partir des FICHIERS produits (pas du code).

Pour chaque sujet x format x étape :
  - fichiers trouvés / attendus (attendu = valeur du YAML, sinon le nombre le
    plus fréquent chez les autres sujets)
  - date du plus ancien et du plus récent fichier
  - ordre : une étape est à jour si son plus ANCIEN fichier est plus récent que
    le plus RÉCENT fichier de chacune de ses dépendances (règle de `make`).
    Le statut "périmé" se propage vers l'aval.
  - les logs SLURM (.out/.err/.sh) ne servent que d'INDICES affichés au survol,
    jamais à décider de la couleur (ils sont absents ou peu fiables pour
    beaucoup de scripts).
  - un historique (JSON) garde la trace des scans successifs.

Usage :
    python pipeline_status.py retinomaps_steps.yaml
    python pipeline_status.py retinomaps_steps.yaml -o status.html --subjects sub-01 sub-02
"""
import argparse
import fnmatch
import glob
import json
import os
import re
import sys
import time
from collections import Counter

try:
    import yaml
except ImportError:
    sys.exit("PyYAML manquant : pip install pyyaml")

HERE = os.path.dirname(os.path.abspath(__file__))
SUSPECT_ERR = re.compile(r"Traceback|CANCELLED|DUE TO TIME LIMIT|oom[-_ ]kill|Out Of Memory|Killed|Segmentation fault", re.I)
END_RE = re.compile(r"End time:\s*([0-9:\- .]+)")
DUR_RE = re.compile(r"Duration:\s*([0-9:\. days,]+)")


# --------------------------------------------------------------------------- config
def load_cfg(path):
    with open(path) as f:
        cfg = yaml.safe_load(f)
    cfg.setdefault("formats", [None])
    cfg.setdefault("tolerance_minutes", 10)
    cfg.setdefault("spread_warning_days", 3)
    cfg.setdefault("ignore", [])
    cfg.setdefault("exclude_subjects", [])
    ids = set()
    for st in cfg["steps"]:
        st.setdefault("label", st["id"])
        st.setdefault("section", "")
        st.setdefault("formats", cfg["formats"])
        st.setdefault("exclude", [])
        deps = []
        for d in st.get("depends", []) or []:
            deps.append({"step": d, "filter": None} if isinstance(d, str) else {"step": d["step"], "filter": d.get("filter")})
        st["depends"] = deps
        ids.add(st["id"])
    for st in cfg["steps"]:
        for d in st["depends"]:
            if d["step"] not in ids:
                sys.exit(f"Étape '{st['id']}' : dépendance inconnue '{d['step']}'")
    return cfg


def as_list(x):
    if x is None:
        return []
    return x if isinstance(x, list) else [x]


def patterns_for(step, fmt):
    out = step["outputs"]
    if isinstance(out, dict):
        out = out.get(fmt) or out.get("default")
    return as_list(out)


# --------------------------------------------------------------------------- scan
_glob_cache = {}


def cached_glob(pattern):
    if pattern not in _glob_cache:
        _glob_cache[pattern] = glob.glob(pattern)
    return _glob_cache[pattern]


def scan_files(root, pats, sub, fmt, ignore, exclude):
    files, unreadable, empty = [], [], []
    seen = set()
    for p in pats:
        full = os.path.join(root, p.format(sub=sub, fmt=fmt or ""))
        for f in cached_glob(full):
            if f in seen:
                continue
            seen.add(f)
            rel = os.path.relpath(f, root)
            if any(fnmatch.fnmatch(rel, ig) for ig in ignore):
                continue
            base = os.path.basename(f)
            if any(fnmatch.fnmatch(base, ex) for ex in exclude):
                continue
            try:
                stt = os.stat(f)
            except OSError:
                unreadable.append(rel)
                continue
            if not os.path.isfile(f):
                continue
            if stt.st_size == 0:
                empty.append(rel)
            files.append((rel, stt.st_mtime, stt.st_size))
    files.sort(key=lambda x: x[1])
    return files, unreadable, empty


def read_tail(path, nbytes=20000):
    try:
        with open(path, "rb") as f:
            f.seek(0, 2)
            size = f.tell()
            f.seek(max(0, size - nbytes))
            return f.read().decode("utf-8", "replace")
    except OSError:
        return ""


def scan_logs(root, pats, sub, fmt):
    """Indices uniquement : dernier job / log, mots-clés suspects, durée si dispo."""
    hits = []
    for p in as_list(pats):
        for f in cached_glob(os.path.join(root, p.format(sub=sub, fmt=fmt or ""))):
            try:
                hits.append((f, os.stat(f).st_mtime))
            except OSError:
                pass
    if not hits:
        return None
    hits.sort(key=lambda x: x[1])
    last_t = hits[-1][1]
    info = {"last": last_t, "n": len(hits), "last_file": os.path.relpath(hits[-1][0], root)}
    # regarde les .out/.err du dernier "lot" (fichiers à moins d'1h du dernier)
    recent = [f for f, t in hits if last_t - t < 3600]
    for f in recent:
        if f.endswith(".out"):
            txt = read_tail(f)
            m, d = END_RE.findall(txt), DUR_RE.findall(txt)
            if m:
                info["end_time"] = m[-1].strip()
            if d:
                info["duration"] = d[-1].strip()
        elif f.endswith(".err"):
            txt = read_tail(f, 6000)
            s = SUSPECT_ERR.findall(txt)
            if s:
                info.setdefault("suspect", []).append(f"{os.path.basename(f)} : {s[-1]}")
    return info


# --------------------------------------------------------------------------- statut
def topo(steps):
    by_id = {s["id"]: s for s in steps}
    order, done, temp = [], set(), set()

    def visit(i):
        if i in done:
            return
        if i in temp:
            sys.exit(f"Cycle de dépendances autour de '{i}'")
        temp.add(i)
        for d in by_id[i]["depends"]:
            visit(d["step"])
        temp.discard(i)
        done.add(i)
        order.append(by_id[i])

    for s in steps:
        visit(s["id"])
    return order


def run(cfg, cfg_path, only_subjects=None, out=None):
    root = cfg["root"]
    tol = cfg["tolerance_minutes"] * 60
    spread = cfg["spread_warning_days"] * 86400
    subs = sorted({os.path.basename(p) for p in glob.glob(os.path.join(root, cfg["subjects_glob"]))
                   if os.path.isdir(p)} - set(cfg["exclude_subjects"]))
    if only_subjects:
        subs = [s for s in subs if s in only_subjects]
    steps = topo(cfg["steps"])
    fmts = cfg["formats"]
    t0 = time.time()

    # 1) scan brut
    raw = {}
    for sub in subs:
        print(f"scan {sub} ...", flush=True)
        for st in steps:
            for fmt in st["formats"]:
                pats = patterns_for(st, fmt)
                if not pats:
                    continue
                files, unread, empty = scan_files(root, pats, sub, fmt, cfg["ignore"], st["exclude"])
                logs = scan_logs(root, st.get("logs"), sub, fmt) if st.get("logs") else None
                raw[(sub, fmt, st["id"])] = dict(files=files, unreadable=unread, empty=empty, logs=logs)

    # 2) nombre attendu (YAML sinon mode chez les sujets)
    expected = {}
    for st in steps:
        for fmt in st["formats"]:
            e = st.get("expected")
            if isinstance(e, dict):
                e = e.get(fmt)
            if e is None:
                counts = [len(raw[(s, fmt, st["id"])]["files"]) for s in subs if (s, fmt, st["id"]) in raw]
                counts = [c for c in counts if c > 0]
                if len(counts) >= 3:
                    e = Counter(counts).most_common(1)[0][0]
                    e = max(e, 0)
                    expected[(fmt, st["id"])] = (e, "inféré")
                    continue
                expected[(fmt, st["id"])] = (None, "")
            else:
                expected[(fmt, st["id"])] = (e, "yaml")

    # 3) statut, dans l'ordre topologique
    res = {}
    for sub in subs:
        for st in steps:
            for fmt in st["formats"]:
                key = (sub, fmt, st["id"])
                if key not in raw:
                    continue
                r = raw[key]
                files = r["files"]
                n = len(files)
                c = dict(n=n, status="missing", flags=[], deps=[])
                exp, exp_src = expected[(fmt, st["id"])]
                c["expected"], c["expected_src"] = exp, exp_src
                if n:
                    c["oldest"], c["newest"] = files[0][1], files[-1][1]
                    c["oldest_file"], c["newest_file"] = files[0][0], files[-1][0]
                # dépendances
                stale, upstream, orphan = [], [], []
                for d in st["depends"]:
                    dfmt = fmt if (sub, fmt, d["step"]) in raw else next(
                        (f for f in fmts if (sub, f, d["step"]) in raw), None)
                    dk = (sub, dfmt, d["step"])
                    if dk not in raw:
                        continue
                    dfiles = raw[dk]["files"]
                    if d["filter"]:
                        dfiles = [x for x in dfiles if fnmatch.fnmatch(os.path.basename(x[0]), d["filter"])]
                    dres = res.get(dk, {})
                    info = dict(step=d["step"], fmt=dfmt, filter=d["filter"], n=len(dfiles),
                                status=dres.get("status", "missing"))
                    if dfiles:
                        info["newest"] = dfiles[-1][1]
                        info["newest_file"] = dfiles[-1][0]
                    c["deps"].append(info)
                    if not n:
                        continue
                    if not dfiles:
                        orphan.append(d["step"])
                    elif dfiles[-1][1] > c["oldest"] + tol:
                        stale.append(d["step"])
                        info["newer"] = True
                    elif dres.get("status") in ("stale", "upstream"):
                        upstream.append(d["step"])
                if n:
                    if stale:
                        c["status"] = "stale"
                    elif upstream:
                        c["status"] = "upstream"
                    elif exp and n < exp:
                        c["status"] = "partial"
                    elif orphan:
                        c["status"] = "orphan"
                    else:
                        c["status"] = "ok"
                    if exp and n > exp:
                        c["flags"].append(f"plus de fichiers que d'habitude ({n} > {exp})")
                    if c["newest"] - c["oldest"] > spread:
                        c["flags"].append("fichiers écrits à des dates très différentes (relance partielle ?)")
                    if orphan:
                        c["flags"].append("dépendance absente : " + ", ".join(orphan))
                if r["unreadable"]:
                    c["flags"].append(f"{len(r['unreadable'])} fichier(s) illisible(s) (permissions)")
                if r["empty"]:
                    c["flags"].append(f"{len(r['empty'])} fichier(s) vide(s)")
                lg = r["logs"]
                if lg:
                    c["logs"] = lg
                    if n and lg["last"] > c["newest"] + max(tol, 1800):
                        c["flags"].append("job/log plus récent que les fichiers (relance en cours ou échouée ?)")
                    if lg.get("suspect"):
                        c["flags"].append("mot-clé suspect dans le .err (à vérifier, pas forcément une erreur)")
                c["files"] = [(f, t) for f, t, _ in files[-300:]]
                c["unreadable"] = r["unreadable"][:50]
                res[key] = c

    # 4) colonnes sans aucun fichier pour aucun sujet -> n/a
    for st in steps:
        for fmt in st["formats"]:
            if all(res.get((s, fmt, st["id"]), {}).get("n", 0) == 0 for s in subs):
                for s in subs:
                    if (s, fmt, st["id"]) in res and not res[(s, fmt, st["id"])]["flags"]:
                        res[(s, fmt, st["id"])]["status"] = "na"

    # 5) historique
    hist_path = os.path.join(os.path.dirname(os.path.abspath(cfg_path)),
                             cfg.get("history_file", "pipeline_status_history.json"))
    try:
        with open(hist_path) as f:
            hist = json.load(f)
    except (OSError, ValueError):
        hist = {}
    now = time.time()
    changes = []
    for (sub, fmt, sid), c in res.items():
        k = f"{sub}|{fmt}|{sid}"
        snap = dict(scan=now, n=c["n"], oldest=c.get("oldest"), newest=c.get("newest"))
        h = hist.setdefault(k, [])
        if not h or (h[-1]["n"], h[-1]["newest"], h[-1]["oldest"]) != (snap["n"], snap["newest"], snap["oldest"]):
            if h:
                changes.append(dict(sub=sub, fmt=fmt, step=sid, before=h[-1], after=snap))
            h.append(snap)
        c["history"] = h[-20:]
    try:
        with open(hist_path, "w") as f:
            json.dump(hist, f)
    except OSError as e:
        print(f"(historique non sauvegardé : {e})")

    # 6) HTML
    data = dict(
        project=cfg.get("project", ""), root=root, scanned=now, duration=time.time() - t0,
        tolerance=cfg["tolerance_minutes"], subjects=subs, formats=[f for f in fmts],
        steps=[dict(id=s["id"], label=s["label"], section=s["section"], script=s.get("script", ""),
                    formats=s["formats"], depends=s["depends"]) for s in cfg["steps"]],
        cells={f"{s}|{f}|{i}": c for (s, f, i), c in res.items()},
        changes=changes[-500:],
    )
    with open(os.path.join(HERE, "template.html")) as f:
        html = f.read().replace("/*__DATA__*/null", json.dumps(data, default=str))
    out = out or os.path.join(os.path.dirname(os.path.abspath(cfg_path)), "pipeline_status.html")
    with open(out, "w") as f:
        f.write(html)
    cnt = Counter(c["status"] for c in res.values())
    print(f"\n{len(subs)} sujets, {time.time()-t0:.0f}s  ->  {out}")
    print("  " + ", ".join(f"{k}: {v}" for k, v in sorted(cnt.items())))
    if changes:
        print(f"  {len(changes)} changement(s) depuis le dernier scan")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("config")
    ap.add_argument("-o", "--output")
    ap.add_argument("--subjects", nargs="*")
    ap.add_argument("--root", help="remplace 'root' du YAML (ex. disque monté)")
    a = ap.parse_args()
    cfg = load_cfg(a.config)
    if a.root:
        cfg["root"] = a.root
    run(cfg, a.config, a.subjects, a.output)
