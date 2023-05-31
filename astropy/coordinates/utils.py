def generate_docs():
    output = []
    try:
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        sc = SkyCoord(1*u.deg, 2*u.deg)
        output =' \n'.join(["sc."+attrname for attrname in dir(sc) if not attrname.startswith('_')])

    except Exception as e:
        output = str(e)

    with open('../../docs/coordinates/TabCompleteList.txt', 'w') as rst_file:
        rst_file.write(output)
# including the line below will make the tox -e build docs command also update TabCompleteList.txt, but runs into persmissions issues on my machine. Unsure if there's a better place to put this
#generate_docs()
