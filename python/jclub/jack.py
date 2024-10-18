from jclub.compiled import jclub_c


def complement0(data):
    return jclub_c.complement0_c(data)


def complement1(data):
    return jclub_c.complement1_c(data)
