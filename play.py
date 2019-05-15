from astropy.io.votable import parse
from astropy.table import Table
import pprint

def parse_votable(filename):
    votable = parse(filename)
    return votable

def show(votobj):
    for resource in votobj.resources:
        for table in resource.tables:
            # ... do something with the table ...
            print(table)  
def use_reader(filename):
    table = Table.read(filename, format='votable')
    return table

if __name__ == '__main__':
    
    from astropy.table import Table
    import pprint
    
    table = Table.read('HLA.169.xml', format='votable')
    pp = pprint.PrettyPrinter(indent=3)
    pp.pprint(table.meta['votable'])
    
    