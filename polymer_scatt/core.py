"""Conservative, deterministic 1D scattering measurements.

All lengths are reciprocal to the input q unit. An accepted Guinier fit is an
apparent scattering size, never proof of an individual polymer-chain size.
Screening thresholds are transparent engineering heuristics, not standards.
"""
from __future__ import annotations

import math
import numpy as np
from scipy.optimize import least_squares
from scipy.signal import find_peaks as _find_peaks, peak_widths, savgol_filter

CRITERIA = {
    "minimum_points": 10, "qmin_rg_max": 0.5, "qmax_rg_max": 1.0,
    "minimum_q_ratio": 1.5, "minimum_guinier_decay": 0.03,
    "maximum_relative_stderr": 0.2, "maximum_endpoint_relative_spread": 0.15,
    "minimum_r2": 0.9, "maximum_reduced_chi2": 3.0,
    "maximum_leading_bins_skipped": 5,
}


def _units(unit):
    unit = str(unit).strip().lower().replace("å", "a").replace("angstrom", "a")
    if unit in {"nm^-1", "nm-1", "1/nm", "nm⁻¹"}:
        return "nm^-1", "nm", 1.0
    if unit in {"a^-1", "a-1", "1/a", "a⁻¹"}:
        return "A^-1", "A", 0.1
    if unit in {"unknown", "", "native"}:
        return "unknown", "inverse native q unit", None
    raise ValueError("q_unit must be unknown, nm^-1 or A^-1; no unit is inferred")


def _linear(x, y, err=None):
    # Scale the abscissa for stable least squares across q unit conversions.
    centre, span = float(np.mean(x)), float(np.ptp(x))
    if span <= 0:
        raise ValueError("The fit has no q span")
    z = (x-centre)/span
    X = np.column_stack([np.ones(len(z)), z])
    w = np.ones(len(z)) if err is None else 1/err
    coef = np.linalg.lstsq(X*w[:, None], y*w, rcond=None)[0]
    residual = y-X@coef
    dof = len(y)-2
    chi = float(np.sum((residual*w)**2))
    cov = np.linalg.inv((X*w[:, None]).T@(X*w[:, None]))
    # Propagate reported sigma as absolute; never shrink it to match residuals.
    cov *= chi/dof if err is None else max(1.0, chi/dof)
    slope = float(coef[1]/span)
    intercept = float(coef[0]-coef[1]*centre/span)
    ss = float(np.sum((y-np.mean(y))**2))
    r2 = 1-float(np.sum(residual**2))/ss if ss > 0 else 0.0
    return slope, intercept, float(math.sqrt(max(cov[1, 1], 0))/span), residual, r2, (chi/dof if err is not None else None)


def _fit_window(q, I, s, ids, start, end, factor):
    qq, ii = q[start:end], I[start:end]
    y, x = np.log(ii), qq**2
    err = None if s is None else s[start:end]/ii
    slope, intercept, slope_se, resid, r2, chi = _linear(x, y, err)
    if slope >= 0:
        return None
    rg = float(math.sqrt(-3*slope))
    rg_se = float(1.5*slope_se/rg)
    qminrg, qmaxrg = float(qq[0]*rg), float(qq[-1]*rg)
    endpoint = []
    trim = max(1, len(qq)//10)
    for lo, hi in [(trim,len(qq)), (0,len(qq)-trim), (trim,len(qq)-trim)]:
        if hi-lo < 6:
            continue
        b = _linear(x[lo:hi],y[lo:hi],None if err is None else err[lo:hi])[0]
        endpoint.append(abs(math.sqrt(-3*b)/rg-1) if b < 0 else 1.0)
    spread = float(max(endpoint, default=1.0))
    # The quadratic coefficient on z in [-.5,.5] measures systematic curvature.
    z = (x-np.mean(x))/np.ptp(x)
    X = np.column_stack([np.ones(len(z)),z,z*z])
    w = np.ones(len(z)) if err is None else 1/err
    beta = np.linalg.lstsq(X*w[:,None],y*w,rcond=None)[0]
    res2 = y-X@beta
    cov2 = np.linalg.pinv((X*w[:,None]).T@(X*w[:,None]))
    cov2 *= float(np.sum((res2*w)**2)/(len(z)-3)) if err is None else max(1.0,float(np.sum((res2*w)**2)/(len(z)-3)))
    se2 = math.sqrt(max(float(cov2[2,2]),0))
    curvature_z = float(abs(beta[2])/max(se2,1e-15))
    curvature_amplitude = float(abs(beta[2])/4)
    decay = float(-slope*np.ptp(x))
    reasons=[]
    if len(qq)<CRITERIA["minimum_points"]: reasons.append("too_few_points")
    if qminrg>CRITERIA["qmin_rg_max"]: reasons.append("insufficient_low_q_access")
    if qmaxrg>CRITERIA["qmax_rg_max"]: reasons.append("outside_guinier_qrg_range")
    if qq[-1]/qq[0]<CRITERIA["minimum_q_ratio"]: reasons.append("narrow_q_interval")
    if decay<CRITERIA["minimum_guinier_decay"]: reasons.append("too_little_resolved_intensity_decay")
    if rg_se/rg>CRITERIA["maximum_relative_stderr"]: reasons.append("poor_size_precision")
    if spread>CRITERIA["maximum_endpoint_relative_spread"]: reasons.append("endpoint_sensitive")
    if r2<CRITERIA["minimum_r2"]: reasons.append("weak_linear_fit")
    if chi is not None and chi>CRITERIA["maximum_reduced_chi2"]: reasons.append("residuals_exceed_reported_errors")
    if curvature_z>3 and curvature_amplitude>max(0.005,0.05*decay): reasons.append("systematic_guinier_curvature")
    if err is not None and np.any(err>0.3): reasons.append("low_signal_log_transform_unreliable")
    return dict(qmin=float(qq[0]),qmax=float(qq[-1]),rg_native=rg,
        rg_nm=None if factor is None else rg*factor,rg_stderr_native=rg_se,
        rg_stderr_nm=None if factor is None else rg_se*factor,qmin_rg=qminrg,
        qmax_rg=qmaxrg,slope=slope,intercept=intercept,slope_stderr=slope_se,
        r2=float(r2),reduced_chi2=chi,n_points=len(qq),start_index=int(ids[start]),
        end_index=int(ids[end-1]),skipped_leading_usable_bins=int(start),
        endpoint_relative_spread=spread,curvature_z=curvature_z,
        curvature_log_amplitude=curvature_amplitude,guinier_log_decay=decay,
        rejection_reasons=reasons)


def analyze_curve(q, intensity, sigma=None, q_unit="unknown", label="", qmin=None, qmax=None):
    """Screen low-q contiguous windows without inventing a Guinier plateau.

    Zero-intensity/zero-sigma pairs are masked detector bins. Other nonpositive
    intensities break log-fit windows; they are never filtered to form a fake
    positive-only curve. Invalid supplied errors are an explicit failure.
    Unknown units allow native-unit Rg only. Indices refer to the original input.
    The fitted background is fixed at zero: subtraction must precede this call.
    """
    unit,length_unit,factor=_units(q_unit)
    q=np.asarray(q,dtype=float).reshape(-1)
    I=np.asarray(intensity,dtype=float).reshape(-1)
    s=None if sigma is None else np.asarray(sigma,dtype=float).reshape(-1)
    out=dict(label=str(label),status="invalid",q_unit=unit,length_unit=length_unit,
        n_input=len(q),n_usable=0,rg_native=None,rg_nm=None,rg_stderr_native=None,
        fit=None,candidate=None,rejection_reasons=[],warnings=[],diagnostics={})
    if unit=="unknown": out["warnings"].append("q units are unconfirmed; lengths cannot be assigned nm or angstrom")
    out["warnings"].append("Apparent scattering Rg is not automatically a polymer-chain Rg; background is fixed at zero")
    if s is None: out["warnings"].append("No measurement errors supplied; fit errors are residual-based and conditional")
    if len(q)!=len(I) or (s is not None and len(s)!=len(q)):
        out["rejection_reasons"]=["array_length_mismatch"]; return out
    if qmin is not None and qmax is not None and qmin>=qmax:
        out["rejection_reasons"]=["invalid_q_bounds"]; return out
    ids=np.arange(len(q))
    masks=(I==0)&(s==0) if s is not None else np.zeros(len(q),dtype=bool)
    keep=~masks
    out["diagnostics"]["zero_zero_masked_bins"]=int(np.sum(masks))
    if qmin is not None: keep &= q>=qmin
    if qmax is not None: keep &= q<=qmax
    q,I,ids=q[keep],I[keep],ids[keep]
    if s is not None: s=s[keep]
    if len(q)<10: out["rejection_reasons"]=["too_few_unmasked_points"]; return out
    if not np.all(np.isfinite(q)) or not np.all(np.isfinite(I)) or np.any(q<=0):
        out["rejection_reasons"]=["nonfinite_values_or_nonpositive_q"]; return out
    if np.any(np.diff(q)<=0): out["rejection_reasons"]=["q_must_be_strictly_increasing"]; return out
    if s is not None and (not np.all(np.isfinite(s)) or np.any(s<=0)):
        out["rejection_reasons"]=["invalid_supplied_sigma"]; return out
    out["n_usable"]=len(q)
    out["diagnostics"].update(qmin_usable=float(q[0]),qmax_usable=float(q[-1]),
        nonpositive_intensity_bins=int(np.sum(I<=0)),background_assumption="zero_after_preprocessing",criteria=CRITERIA.copy())
    candidates=[]
    # Fixed, short search: do not roam to an unrelated high-q shoulder.
    for start in (0,1,2,3,5):
        if start>=len(q) or I[start]<=0: continue
        bad=np.flatnonzero((I[start:]<=0))
        stop=start+int(bad[0]) if len(bad) else len(q)
        # Mask gaps inside a curve also break a contiguous fit region.
        gaps=np.flatnonzero(np.diff(ids[start:stop])>1)
        if len(gaps): stop=min(stop,start+int(gaps[0])+1)
        available=stop-start
        if available<10: continue
        # Log-spaced candidate lengths cover the entire positive low-q segment.
        lengths=np.unique(np.rint(np.geomspace(10,available,min(12,available-9))).astype(int))
        for n in lengths:
            try: f=_fit_window(q,I,s,ids,start,start+int(n),factor)
            except (ValueError,np.linalg.LinAlgError,FloatingPointError): continue
            if f is not None: candidates.append(f)
    out["diagnostics"]["candidate_windows_evaluated"]=len(candidates)
    if not candidates:
        out["status"]="rejected";out["rejection_reasons"]=["no_contiguous_positive_low_q_window_with_negative_slope"];return out
    def rank(f):
        # Prefer gate-valid fits, then broad information with early low-q access.
        penalty=len(f["rejection_reasons"])
        continuous=max(f["qmin_rg"]/0.5-1,0)+max(f["qmax_rg"]-1,0)
        information=f["guinier_log_decay"]/(1+5*f["endpoint_relative_spread"])
        return (penalty, continuous, -information, f["skipped_leading_usable_bins"])
    best=min(candidates,key=rank)
    out["candidate"]=best
    out["rejection_reasons"]=best["rejection_reasons"].copy()
    if best["skipped_leading_usable_bins"]:
        out["warnings"].append(f"Fit skips {best['skipped_leading_usable_bins']} leading usable bins; inspect their residuals")
    if not best["rejection_reasons"]:
        out.update(status="accepted",fit=best,rg_native=best["rg_native"],rg_nm=best["rg_nm"],rg_stderr_native=best["rg_stderr_native"])
    else: out["status"]="rejected"
    return out


def debye_form_factor(q, rg):
    """Stable monodisperse Gaussian-coil form factor, normalized at q=0."""
    z=(np.asarray(q,dtype=float)*rg)**2
    small=z<1e-3
    answer=np.empty_like(z)
    answer[small]=1-z[small]/3+z[small]**2/12-z[small]**3/60
    answer[~small]=2*(np.expm1(-z[~small])+z[~small])/z[~small]**2
    return answer


def fit_debye(q,intensity,sigma=None,q_unit="unknown",qmin=None,qmax=None,fit_background=False):
    """Optional Gaussian-chain hypothesis, not a validated chain-size result.

    Background is fixed to zero unless explicitly requested. Negative measured
    intensities are retained in this intensity-space fit.
    """
    unit,length_unit,factor=_units(q_unit)
    q=np.asarray(q,dtype=float);I=np.asarray(intensity,dtype=float)
    s=None if sigma is None else np.asarray(sigma,dtype=float)
    if q.shape!=I.shape or (s is not None and q.shape!=s.shape): raise ValueError("Array lengths differ")
    mask=~((I==0)&(s==0)) if s is not None else np.ones(q.shape,dtype=bool)
    if qmin is not None:mask&=q>=qmin
    if qmax is not None:mask&=q<=qmax
    q,I=q[mask],I[mask]
    if s is not None:s=s[mask]
    if len(q)<10 or not np.all(np.isfinite(q)) or not np.all(np.isfinite(I)) or np.any(q<=0) or np.any(np.diff(q)<=0):raise ValueError("Need >=10 finite points on strictly increasing positive q")
    if s is not None and (np.any(s<=0) or not np.all(np.isfinite(s))):raise ValueError("Invalid supplied sigma")
    scale=max(float(np.max(np.abs(I))),1e-30)
    Y=I/scale;err=np.ones(len(q)) if s is None else s/scale
    lower=[-30,math.log(1e-4/q[-1])];upper=[30,math.log(1e4/q[0])]
    x0=[math.log(max(float(np.max(Y)),1e-5)), math.log(1/np.sqrt(q[0]*q[-1]))]
    if fit_background: lower+=[-10];upper+=[10];x0+=[0]
    def model(p):return math.exp(p[0])*debye_form_factor(q,math.exp(p[1]))+(p[2] if fit_background else 0)
    fit=least_squares(lambda p:(model(p)-Y)/err,x0,bounds=(lower,upper),max_nfev=500)
    rg=math.exp(float(fit.x[1]));dof=len(q)-len(fit.x)
    rss=float(np.sum(fit.fun**2));normal=fit.jac.T@fit.jac
    condition=float(np.linalg.cond(normal));cov=np.linalg.pinv(normal)
    cov*=rss/dof if s is None else max(1,rss/dof)
    se=rg*math.sqrt(max(float(cov[1,1]),0))
    corr=float(cov[0,1]/math.sqrt(max(float(cov[0,0]*cov[1,1]),1e-300)))
    reasons=[]
    if not fit.success:reasons.append("optimizer_did_not_converge")
    if q[0]*rg>0.5:reasons.append("low_q_plateau_not_resolved")
    if q[-1]*rg<1:reasons.append("coil_crossover_not_resolved")
    if se/rg>0.2:reasons.append("poor_parameter_precision")
    if abs(corr)>0.995:reasons.append("amplitude_size_degeneracy")
    if not np.isfinite(condition) or condition>1e12:reasons.append("ill_conditioned_covariance")
    if s is not None and rss/dof>3:reasons.append("model_mismatch_to_reported_errors")
    if any(abs(fit.x[i]-lower[i])<1e-3 or abs(fit.x[i]-upper[i])<1e-3 for i in range(len(fit.x))):reasons.append("parameter_at_bound")
    return dict(status="hypothesis_fit" if not reasons else "unresolved",model="monodisperse_Debye_Gaussian_coil",q_unit=unit,length_unit=length_unit,
        rg_native=rg,rg_nm=None if factor is None else rg*factor,rg_stderr_native=se,
        amplitude=float(math.exp(fit.x[0])*scale),background=float(fit.x[2]*scale) if fit_background else 0.0,
        qmin=float(q[0]),qmax=float(q[-1]),qmin_rg=float(q[0]*rg),qmax_rg=float(q[-1]*rg),
        reduced_chi2=None if s is None else rss/dof,amplitude_rg_correlation=corr,
        covariance_condition_number=condition if np.isfinite(condition) else None,n_points=len(q),
        rejection_reasons=reasons,warning="Model-dependent apparent coil size; metadata and contrast must support individual Gaussian-chain scattering")


def find_peaks(q,intensity,q_unit="unknown",prominence_fraction=0.05,smooth_points=11,sigma=None):
    """Descriptive peak positions/widths; no phase or crystallinity inference.

    If sigma is supplied, zero-intensity/zero-sigma pairs are masked. Contiguous
    valid runs are processed separately so smoothing cannot bridge a detector
    gap. Without sigma, measured zeros are retained, not assumed to be masks.
    Errors validate masks only; this peak detector is not uncertainty-weighted.
    """
    unit,length_unit,factor=_units(q_unit)
    q=np.asarray(q,dtype=float);I=np.asarray(intensity,dtype=float)
    if q.shape!=I.shape or len(q)<5 or not np.all(np.isfinite(q)) or not np.all(np.isfinite(I)) or np.any(q<=0) or np.any(np.diff(q)<=0):
        raise ValueError("Need finite intensities on strictly increasing positive q")
    if not 0<prominence_fraction<1:raise ValueError("prominence_fraction must lie between zero and one")
    if int(smooth_points)<1:raise ValueError("smooth_points must be positive")
    s=None if sigma is None else np.asarray(sigma,dtype=float)
    if s is not None and s.shape!=q.shape:raise ValueError("Sigma shape differs from q")
    mask=(I==0)&(s==0) if s is not None else np.zeros(len(q),dtype=bool)
    if s is not None and (not np.all(np.isfinite(s[~mask])) or np.any(s[~mask]<=0)):
        raise ValueError("Invalid supplied sigma")
    ids=np.flatnonzero(~mask)
    result=[];windows=[]
    for indices in np.split(ids,np.flatnonzero(np.diff(ids)>1)+1):
        if len(indices)<5:continue
        qq,yy=q[indices],I[indices]
        n=min(int(smooth_points)//2*2+1,len(qq)//2*2-1)
        windows.append(n)
        sm=savgol_filter(yy,n,2) if n>=5 else yy.copy()
        # Suppress floating-point Savitzky-Golay ripples on a flat baseline.
        prominence=max(float(np.ptp(sm))*prominence_fraction,
            float(np.max(np.abs(sm)))*np.finfo(float).eps*100,1e-30)
        peaks,props=_find_peaks(sm,prominence=prominence)
        widths=peak_widths(sm,peaks,rel_height=.5) if len(peaks) else None
        for k,local in enumerate(peaks):
            idx=indices[local]
            left=float(np.interp(widths[2][k],np.arange(len(qq)),qq));right=float(np.interp(widths[3][k],np.arange(len(qq)),qq))
            d=2*math.pi/float(q[idx])
            result.append(dict(index=int(idx),q_peak=float(q[idx]),d_native=d,d_nm=None if factor is None else d*factor,
                intensity=float(I[idx]),prominence=float(props["prominences"][k]),half_prominence_width_q=right-left))
    return dict(q_unit=unit,length_unit=length_unit,peaks=result,
        smoothing_points=max(windows,default=0),smoothing_points_by_segment=windows,
        zero_zero_masked_bins=int(np.sum(mask)),
        warning="Descriptive peak spacings only; widths are half-prominence widths, not instrument-corrected crystallite sizes. Smoothing uses point indices; irregular q grids require review.")


peak_analysis=find_peaks
