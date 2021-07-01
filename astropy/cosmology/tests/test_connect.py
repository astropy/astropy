# # Licensed under a 3-clause BSD style license - see LICENSE.rst
#
# import os
#
# import pytest
#
# from astropy import cosmology
# from astropy.cosmology import Cosmology
# from astropy.table import QTable, vstack
# from astropy.cosmology.connect import CosmologyRead
# from astropy.utils.exceptions import AstropyUserWarning
#
# cosmo_instances = cosmology.parameters.available
# save_formats = ["json", "ascii.ecsv"]
#
#
# # make a common directory for reading / writing cosmologies
# @pytest.fixture(scope="session")
# def cosmo_dir(tmpdir_factory):
#     drct = tmpdir_factory.mktemp("cosmo")
#     return drct
#
#
# # -----------------------------------------------------------------------------
#
#
# class TestReadWriteCosmology:
#
#     def test_instantiate_read(self):
#         # no error on base class
#         assert isinstance(Cosmology.read, CosmologyRead)
#
#         # Warns when initialized. Cannot be used.
#         with pytest.warns(AstropyUserWarning):
#             assert cosmology.realizations.Planck18.read is NotImplemented
#
#     @pytest.mark.parametrize("format", save_formats)
#     @pytest.mark.parametrize("instance", cosmo_instances)
#     def test_write_then_read_file(self, cosmo_dir, instance, format):
#         """Read tests happen later."""
#         cosmo = getattr(cosmology.realizations, instance)
#         fname = cosmo_dir / f"{instance}.{format}"
#
#         cosmo.write(str(fname), format=format)
#
#         # Also test kwarg "overwrite"
#         assert os.path.exists(str(fname))  # file exists
#         with pytest.raises(IOError):
#             cosmo.write(str(fname), format=format, overwrite=False)
#
#         assert os.path.exists(str(fname))  # overwrite file existing file
#         cosmo.write(str(fname), format=format, overwrite=True)
#
#         # Read back
#         got = Cosmology.read(cosmo_dir / f"{instance}.{format}", format=format)
#
#         assert got.name == cosmo.name
#         assert got == cosmo
#
#     @pytest.mark.parametrize("instance", cosmo_instances)
#     def test_to_mapping_instance(self, instance):
#         instance = getattr(cosmology.realizations, instance)
#         m = instance.write.to_mapping()
#
#         assert isinstance(m, dict)
#         assert "cosmology" in m
#
#     @pytest.mark.parametrize("instance", cosmo_instances)
#     def test_to_table_instance(self, instance):
#         instance = getattr(cosmology.realizations, instance)
#         t = instance.write.to_table()
#
#         assert isinstance(t, QTable)
#         assert "cosmology" in t.meta
#
#     @pytest.mark.parametrize("instance", cosmo_instances)
#     def test_from_mapping_instance(self, instance):
#         expected = getattr(cosmology.realizations, instance)
#         params = getattr(cosmology.parameters, instance)
#
#         # read with mismatching parameters errors
#         with pytest.raises(TypeError, match="there are unused parameters"):
#             Cosmology.read.from_mapping(params)
#
#         # unless mismatched are moved to meta
#         got = Cosmology.read.from_mapping(params, move_to_meta=True)
#         assert got == expected
#
#     @pytest.mark.parametrize("format", save_formats[1:])  # skip json
#     @pytest.mark.parametrize("instance", cosmo_instances)
#     def test_from_table_instance(self, cosmo_dir, instance, format):
#         expected = getattr(cosmology.realizations, instance)
#
#         tbl = QTable.read(cosmo_dir / f"{instance}.{format}", format=format)
#         got = Cosmology.read.from_table(tbl)
#         assert got == expected
#
#
# # -----------------------------------------------------------------------------
#
#
# @pytest.mark.parametrize("instance", cosmo_instances)
# def test_round_trip_of_mapping_instance(instance):
#     expected = getattr(cosmology.realizations, instance)
#
#     got = Cosmology.read.from_mapping(expected.write.to_mapping())
#     assert got == expected
#
#
# @pytest.mark.parametrize("instance", cosmo_instances)
# class Test_round_trip_of_table_instance:
#
#     def test_1D(self, instance):
#         expected = getattr(cosmology.realizations, instance)
#
#         got = Cosmology.read.from_table(expected.write.to_table())
#         assert got == expected
#
#     def test_multirow(self, instance):
#         expected = getattr(cosmology.realizations, instance)
#         expected1 = expected.clone(name="Other")
#         t = vstack([expected.write.to_table(), expected1.write.to_table()])
#
#         # error for no index
#         with pytest.raises(ValueError, match="multi-row"):
#             Cosmology.read.from_table(t)
#
#         got = Cosmology.read.from_table(t, index=0)
#         assert got == expected
#
#         got = Cosmology.read.from_table(t, index=1)
#         assert got == expected1
#
#         # Now test index is string
#         got = Cosmology.read.from_table(t, index="Other")
#         assert got == expected1
