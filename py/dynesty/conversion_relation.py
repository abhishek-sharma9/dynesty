import numpy as np
import lal

def _a0(f_lower):
    """Used in calculating chirp times: see Cokelaer, arxiv.org:0706.4437
       appendix 1, also lalinspiral/python/sbank/tau0tau3.py.
    """
    return 5. / (256. * (np.pi * f_lower)**(8./3.))


def _a3(f_lower):
    """Another parameter used for chirp times"""
    return np.pi / (8. * (np.pi * f_lower)**(5./3.))

def mtotal_eta_from_tau0_tau3(tau0, tau3, f_lower,
                          in_seconds=False):
    r"""Returns total mass from :math:`\tau_0, \tau_3`."""
    mtotal = (tau3 / _a3(f_lower)) / (tau0 / _a0(f_lower))
    eta = mtotal**(-2./3.) * (_a3(f_lower) / tau3)
    if not in_seconds:
        # convert back to solar mass units
        mtotal /= lal.MTSUN_SI
    return mtotal, eta

def mass1_mass2_from_mtotal_eta(m_total, eta):
    mass1 = 0.5 * m_total * (1.0 + (1.0 - 4.0 * eta)**0.5)
    mass2 = 0.5 * m_total * (1.0 - (1.0 - 4.0 * eta)**0.5)
    return mass1, mass2
