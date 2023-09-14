"""
Created on Feb 26, 2021

@author: laurentmichel
"""
import os

from astropy.io.votable import parse
from astropy.io.votable.tree import MivotBlock, Resource, VOTableFile

data_path = os.path.dirname(os.path.realpath(__file__))

"""
Read out Test
"""

# Read a valid VOTable
# No invalid element detected,
# The block is returned as an XML String
vpath = os.path.join(data_path, "test.7.xml")
votable = parse(vpath)
for resource in votable.resources:
    print(f"Resource type: {resource.type}  Mivot Block: {resource.mivot_block}")

    for resource_meta in resource.resources:
        print(f"Resource type: {resource_meta.type} Mivot Block:")
        print(resource_meta.mivot_block)

# Read an invalid VOTable
# An unexpected element has been found in the mivot block
# The mivot block pointer returns a block with just REPORT in error
vpath = os.path.join(data_path, "test.7.ko.xml")
votable = parse(vpath)
for resource in votable.resources:
    print(f"Resource type: {resource.type}  Mivot Block: {resource.mivot_block}")

    for resource_meta in resource.resources:
        print(f"Resource type: {resource_meta.type} Mivot Block:")
        print(resource_meta.mivot_block)

"""
Write Test
"""
vpath = os.path.join(data_path, "test.7.out.xml")

# Create am empty VOTable
votable = VOTableFile()
# Create the resource that will host both data table and mivot resource.
resource = Resource()
resource.type = "results"
# Create the resource that will host the mivot.
meta_resource = Resource()
meta_resource.type = "meta"
# A dummy mivot block for the test.
resource.resources.append(meta_resource)
mivot_block = MivotBlock("""
<dm-mapping:VODML xmlns:dm-mapping="http://www.ivoa.net/xml/merged-syntax" >
  <dm-mapping:REPORT/>
  <dm-mapping:GLOBALS/>
</dm-mapping:VODML>
"""
                         )
# Add the mivot resource
meta_resource.mivot_block = mivot_block
votable.resources.append(resource)
# Save the VOTable
votable.to_xml(vpath)
# and read it again to retrieve the mivot
with open(vpath) as result:
    print(result.read())
