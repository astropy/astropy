# Type stub for scalar_inv_efuncs.pyx. Keep in sync with its signatures.

def lcdm_inv_efunc_norel(z: float, Om0: float, Ode0: float, Ok0: float) -> float: ...
def lcdm_inv_efunc_nomnu(
    z: float,
    Om0: float,
    Ode0: float,
    Ok0: float,
    Or0: float,
) -> float: ...
def lcdm_inv_efunc(
    z: float,
    Om0: float,
    Ode0: float,
    Ok0: float,
    Ogamma0: float,
    NeffPerNu: float,
    nmasslessnu: int,
    nu_y: list[float],
) -> float: ...
def flcdm_inv_efunc_norel(z: float, Om0: float, Ode0: float) -> float: ...
def flcdm_inv_efunc_nomnu(z: float, Om0: float, Ode0: float, Or0: float) -> float: ...
def flcdm_inv_efunc(
    z: float,
    Om0: float,
    Ode0: float,
    Ogamma0: float,
    NeffPerNu: float,
    nmasslessnu: int,
    nu_y: list[float],
) -> float: ...
def wcdm_inv_efunc_norel(
    z: float,
    Om0: float,
    Ode0: float,
    Ok0: float,
    w0: float,
) -> float: ...
def wcdm_inv_efunc_nomnu(
    z: float,
    Om0: float,
    Ode0: float,
    Ok0: float,
    Or0: float,
    w0: float,
) -> float: ...
def wcdm_inv_efunc(
    z: float,
    Om0: float,
    Ode0: float,
    Ok0: float,
    Ogamma0: float,
    NeffPerNu: float,
    nmasslessnu: int,
    nu_y: list[float],
    w0: float,
) -> float: ...
def fwcdm_inv_efunc_norel(z: float, Om0: float, Ode0: float, w0: float) -> float: ...
def fwcdm_inv_efunc_nomnu(
    z: float,
    Om0: float,
    Ode0: float,
    Or0: float,
    w0: float,
) -> float: ...
def fwcdm_inv_efunc(
    z: float,
    Om0: float,
    Ode0: float,
    Ogamma0: float,
    NeffPerNu: float,
    nmasslessnu: int,
    nu_y: list[float],
    w0: float,
) -> float: ...
def w0wacdm_inv_efunc_norel(
    z: float,
    Om0: float,
    Ode0: float,
    Ok0: float,
    w0: float,
    wa: float,
) -> float: ...
def w0wacdm_inv_efunc_nomnu(
    z: float,
    Om0: float,
    Ode0: float,
    Ok0: float,
    Or0: float,
    w0: float,
    wa: float,
) -> float: ...
def w0wacdm_inv_efunc(
    z: float,
    Om0: float,
    Ode0: float,
    Ok0: float,
    Ogamma0: float,
    NeffPerNu: float,
    nmasslessnu: int,
    nu_y: list[float],
    w0: float,
    wa: float,
) -> float: ...
def fw0wacdm_inv_efunc_norel(
    z: float,
    Om0: float,
    Ode0: float,
    w0: float,
    wa: float,
) -> float: ...
def fw0wacdm_inv_efunc_nomnu(
    z: float,
    Om0: float,
    Ode0: float,
    Or0: float,
    w0: float,
    wa: float,
) -> float: ...
def fw0wacdm_inv_efunc(
    z: float,
    Om0: float,
    Ode0: float,
    Ogamma0: float,
    NeffPerNu: float,
    nmasslessnu: int,
    nu_y: list[float],
    w0: float,
    wa: float,
) -> float: ...
def wpwacdm_inv_efunc_norel(
    z: float,
    Om0: float,
    Ode0: float,
    Ok0: float,
    wp: float,
    apiv: float,
    wa: float,
) -> float: ...
def wpwacdm_inv_efunc_nomnu(
    z: float,
    Om0: float,
    Ode0: float,
    Ok0: float,
    Or0: float,
    wp: float,
    apiv: float,
    wa: float,
) -> float: ...
def wpwacdm_inv_efunc(
    z: float,
    Om0: float,
    Ode0: float,
    Ok0: float,
    Ogamma0: float,
    NeffPerNu: float,
    nmasslessnu: int,
    nu_y: list[float],
    wp: float,
    apiv: float,
    wa: float,
) -> float: ...
def fwpwacdm_inv_efunc_norel(
    z: float,
    Om0: float,
    Ode0: float,
    wp: float,
    apiv: float,
    wa: float,
) -> float: ...
def fwpwacdm_inv_efunc_nomnu(
    z: float,
    Om0: float,
    Ode0: float,
    Or0: float,
    wp: float,
    apiv: float,
    wa: float,
) -> float: ...
def fwpwacdm_inv_efunc(
    z: float,
    Om0: float,
    Ode0: float,
    Ogamma0: float,
    NeffPerNu: float,
    nmasslessnu: int,
    nu_y: list[float],
    wp: float,
    apiv: float,
    wa: float,
) -> float: ...
def w0wzcdm_inv_efunc_norel(
    z: float,
    Om0: float,
    Ode0: float,
    Ok0: float,
    w0: float,
    wz: float,
) -> float: ...
def w0wzcdm_inv_efunc_nomnu(
    z: float,
    Om0: float,
    Ode0: float,
    Ok0: float,
    Or0: float,
    w0: float,
    wz: float,
) -> float: ...
def w0wzcdm_inv_efunc(
    z: float,
    Om0: float,
    Ode0: float,
    Ok0: float,
    Ogamma0: float,
    NeffPerNu: float,
    nmasslessnu: int,
    nu_y: list[float],
    w0: float,
    wz: float,
) -> float: ...
def fw0wzcdm_inv_efunc_norel(
    z: float,
    Om0: float,
    Ode0: float,
    w0: float,
    wz: float,
) -> float: ...
def fw0wzcdm_inv_efunc_nomnu(
    z: float,
    Om0: float,
    Ode0: float,
    Or0: float,
    w0: float,
    wz: float,
) -> float: ...
def fw0wzcdm_inv_efunc(
    z: float,
    Om0: float,
    Ode0: float,
    Ogamma0: float,
    NeffPerNu: float,
    nmasslessnu: int,
    nu_y: list[float],
    w0: float,
    wz: float,
) -> float: ...
