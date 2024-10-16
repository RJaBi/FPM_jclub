from .compiled import jclub_c


def plaquette(data, mustart, muend, nuend):
    return jclub_c.plaquette_c(data, mustart, muend, nuend)


def plaquette_kolaorder(data, mustart, muend, nuend):
    return jclub_c.plaquette_kolaorder_c(data, mustart, muend, nuend)


def polyakov(data):
    return jclub_c.polyakov_c(data)


def polyakov_kolaorder(data):
    return jclub_c.poly_kolaorder_c(data)
