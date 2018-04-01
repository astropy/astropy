# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import numpy as np
from numpy.testing import assert_array_equal

from asdf import yamlutil

from astropy import table
from astropy.io import fits
from ...types import AstropyAsdfType


class FitsType(AstropyAsdfType):
    name = 'fits/fits'
    types = ['astropy.io.fits.HDUList']
    requires = ['astropy']

    @classmethod
    def from_tree(cls, data, ctx):
        hdus = []
        first = True
        for hdu_entry in data:
            header = fits.Header([fits.Card(*x) for x in hdu_entry['header']])
            data = hdu_entry.get('data')
            if data is not None:
                try:
                    data = data.__array__()
                except ValueError:
                    data = None
            if first:
                hdu = fits.PrimaryHDU(data=data, header=header)
                first = False
            elif data.dtype.names is not None:
                hdu = fits.BinTableHDU(data=data, header=header)
            else:
                hdu = fits.ImageHDU(data=data, header=header)
            hdus.append(hdu)
        hdulist = fits.HDUList(hdus)
        return hdulist

    @classmethod
    def to_tree(cls, hdulist, ctx):
        units = []
        for hdu in hdulist:
            header_list = []
            for card in hdu.header.cards:
                if card.comment:
                    new_card = [card.keyword, card.value, card.comment]
                else:
                    if card.value:
                        new_card = [card.keyword, card.value]
                    else:
                        if card.keyword:
                            new_card = [card.keyword]
                        else:
                            new_card = []
                header_list.append(new_card)

            hdu_dict = {}
            hdu_dict['header'] = header_list
            if hdu.data is not None:
                if hdu.data.dtype.names is not None:
                    data = table.Table(hdu.data)
                else:
                    data = hdu.data
                hdu_dict['data'] = yamlutil.custom_tree_to_tagged_tree(data, ctx)

            units.append(hdu_dict)

        return units

    @classmethod
    def reserve_blocks(cls, data, ctx):
        for hdu in data:
            if hdu.data is not None:
                yield ctx.blocks.find_or_create_block_for_array(hdu.data, ctx)

    @classmethod
    def assert_equal(cls, old, new):
        for hdua, hdub in zip(old, new):
            assert_array_equal(hdua.data, hdub.data)
            for carda, cardb in zip(hdua.header.cards, hdub.header.cards):
                assert tuple(carda) == tuple(cardb)
