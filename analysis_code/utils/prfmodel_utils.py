import numpy as np
from keras import ops
from prfmodel.models.base import BaseTuning
from prfmodel.models.prf import predict_gaussian_response
from prfmodel.stimuli import PRFStimulus, PRFStimulusTensors
from prfmodel.fitters.adapter import ParameterTransform
 
 
class PolarGaussian2DPRFTuning(BaseTuning[PRFStimulus, PRFStimulusTensors]):
    """
    2D isotropic Gaussian pRF tuning model parameterized in polar coordinates.
 
    Identical to prfmodel's Gaussian2DPRFTuning, except that the pRF center is
    given as eccentricity and polar angle instead of x and y. Used for the grid
    search, so that GridFitter's cartesian product over (ecc, polar, sigma)
    reproduces the prfpy polar grid.
 
    Parameters (columns of the parameter DataFrame)
    ----------
    ecc : pRF center eccentricity (dva)
    polar : pRF center polar angle (rad), 0 = right horizontal meridian
    sigma : pRF size (dva)
    """
 
    @property
    def parameter_names(self):
        return ["ecc", "polar", "sigma"]
 
    def call(self, stimulus, parameters):
        ecc = parameters["ecc"]
        polar = parameters["polar"]
        # y first, x second: same order as the stimulus grid
        mu = ops.stack([ecc * ops.sin(polar), ecc * ops.cos(polar)], axis=1)
        return predict_gaussian_response(stimulus.grid, mu, parameters[["sigma"]])
 
 
def polar_to_cartesian(params):
    """
    Convert pRF centers from (ecc, polar) to (mu_x, mu_y).
 
    Parameters
    ----------
    params : pandas.DataFrame
        Parameters with columns 'ecc' and 'polar' (e.g. output of a grid fit
        with PolarGaussian2DPRFTuning).
 
    Returns
    -------
    params_cart : pandas.DataFrame
        Same parameters with 'mu_x' and 'mu_y' instead of 'ecc' and 'polar'.
    """
    params_cart = params.assign(mu_x=params["ecc"] * np.cos(params["polar"]),
                                mu_y=params["ecc"] * np.sin(params["polar"]))
    return params_cart.drop(columns=["ecc", "polar"])
 

def bounded_transform(parameter_names, lower, upper):
    """
    Keep parameters inside (lower, upper) during SGD fitting.
 
    The optimizer works on an unbounded scale (logit); before every prediction
    the value is mapped back with a scaled sigmoid, so the model only ever sees
    values inside the bounds (equivalent to prfpy's bounded iterative fit).
 
    Parameters
    ----------
    parameter_names : list of str
        Parameters that share these bounds.
    lower, upper : float
        Open bounds. Start values must lie strictly inside (see clip_to_bounds).
 
    Returns
    -------
    transform : prfmodel.fitters.adapter.ParameterTransform
        Transform to put inside an Adapter.
    """
    return ParameterTransform(
        parameter_names,
        lambda p: ops.log((p - lower) / (upper - p)),          # natural -> unbounded
        lambda u: lower + (upper - lower) * ops.sigmoid(u))    # unbounded -> natural
 
 
def clip_to_bounds(params, bounds, margin=1e-3):
    """
    Clip start values strictly inside open bounds (required by bounded_transform).
 
    Parameters
    ----------
    params : pandas.DataFrame
        Start parameters.
    bounds : dict
        {parameter_name: (lower, upper)}.
    margin : float, optional
        Distance from the bounds as a fraction of the bound range.
 
    Returns
    -------
    params_clipped : pandas.DataFrame
        Copy of params with the bounded columns clipped.
    """
    params_clipped = params.copy()
    for name, (lower, upper) in bounds.items():
        eps = margin * (upper - lower)
        params_clipped[name] = params_clipped[name].clip(lower + eps, upper - eps)
    return params_clipped
 