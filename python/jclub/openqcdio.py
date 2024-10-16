from jclub.compiled import jclub_c


def readopenqcd(filename, nx, ny, nz, nt):
    return jclub_c.readgaugefieldunopenqcd_c(filename, nx, ny, nz, nt)


def readopenqcd_kolaorder(filename, nx, ny, nz, nt):
    return jclub_c.readgaugefieldunopenqcdunkolaorder_c(filename, nx, ny, nz, nt)


#def say_hello():
#    {{cookiecutter.project_name}}_c.say_hello_c()
