

import numpy as np
import pytest
from astropy.utils.metaclasses import FactoryMeta


class FactoryMeta_TestBase:
    pass


class Test_SupportsAsType(FactoryMeta_TestBase):

    @pytest.fixture
    def inmixbase(self):

        class SupportsAsType(metaclass=FactoryMeta):

            @classmethod
            def _inmix_make_instance(cls, data, /, *args, **kwargs):
                inmixcls = cls._inmix_get_subclass(type(data))
                return inmixcls(data, *args, **kwargs)

            def astype(self, dtype):
                return type(self)([dtype(x) for x in self])

        return SupportsAsType

    # ===============================================================

    def test_defined_cls(self, inmixbase):

        class SupportsAsTypeList(list, inmixbase, data_cls=list):
            pass

        assert issubclass(SupportsAsTypeList, list)
        assert issubclass(SupportsAsTypeList, inmixbase)

        inst = SupportsAsTypeList([1, 2, 3])
        assert isinstance(inst, list)
        assert isinstance(inst, inmixbase)

        # test values
        assert np.array_equal(inst, [1, 2, 3])
        assert all(isinstance(x, int) for x in inst)

        # test added method
        assert hasattr(inst, "astype")
        assert all(isinstance(x, float) for x in inst.astype(float))

    @pytest.mark.parametrize("kind", [list, tuple])
    def test_class_factory(self, inmixbase, kind):

        cls1 = inmixbase(kind)
        cls2 = inmixbase(kind)

        assert cls1 is cls2

        assert issubclass(cls1, kind)
        assert issubclass(cls1, inmixbase)

        inst = cls1([1, 2, 3])
        assert isinstance(inst, kind)
        assert isinstance(inst, inmixbase)

        # test values
        assert np.array_equal(inst, [1, 2, 3])
        assert all(isinstance(x, int) for x in inst)

        # test added method
        assert hasattr(inst, "astype")
        assert all(isinstance(x, float) for x in inst.astype(float))

    @pytest.mark.parametrize("kind", [list, tuple])
    def test_instance_factory(self, inmixbase, kind):

        obj = kind([1, 2, 3])
        inst = inmixbase(obj)
        assert isinstance(inst, kind)
        assert isinstance(inst, inmixbase)

        # test values
        assert np.array_equal(inst, obj)
        assert all(isinstance(x, int) for x in inst)

        # test added method
        assert hasattr(inst, "astype")
        assert all(isinstance(x, float) for x in inst.astype(float))
