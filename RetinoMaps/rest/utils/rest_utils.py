def compute_winners(df_corr: pd.DataFrame) -> np.ndarray:
    """
    Return a 1-D array of winning seed labels (1-based int) per parcel
    Self-seed parcels are masked to NaN before the argmax
    Parcels where all seeds are NaN receive NaN
    """
    df = df_corr.copy()

    for seed, plist in seed_to_parcels.items():
        for p in plist:
            if seed in df.index and p in df.columns:
                df.loc[seed, p] = np.nan

    winners = []
    for parcel in df.columns:
        col = df[parcel]
        winners.append(np.nan if col.isna().all() else seed_to_number[col.idxmax()])

    return np.array(winners)
