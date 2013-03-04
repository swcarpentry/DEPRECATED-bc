"""
ipythonblocks provides a BlockGrid class that displays a colored grid in the
IPython Notebook. The colors can be manipulated, making it useful for
practicing control flow stuctures and quickly seeing the results.

"""

# This file is copyright 2013 by Matt Davis and covered by the license at
# https://github.com/jiffyclub/ipythonblocks/blob/master/LICENSE.txt

import copy
import itertools
import numbers
import os
import sys
import time
import uuid

from operator import iadd

from IPython.display import HTML, display, clear_output

if sys.version_info[0] >= 3:
    xrange = range
    from functools import reduce

__all__ = ('Block', 'BlockGrid', 'Pixel', 'ImageGrid',
           'InvalidColorSpec', 'show_color', 'embed_colorpicker',
           'colors', '__version__')
__version__ = '1.5'

_TABLE = ('<style type="text/css">'
          'table.blockgrid {{border: none;}}'
          ' .blockgrid tr {{border: none;}}'
          ' .blockgrid td {{padding: 0px;}}'
          ' #blocks{0} td {{border: {1}px solid white;}}'
          '</style>'
          '<table id="blocks{0}" class="blockgrid"><tbody>{2}</tbody></table>')
_TR = '<tr>{0}</tr>'
_TD = ('<td title="{0}" style="width: {1}px; height: {1}px;'
       'background-color: {2};"></td>')
_RGB = 'rgb({0}, {1}, {2})'
_TITLE = 'Index: [{0}, {1}]&#10;Color: ({2}, {3}, {4})'

_SINGLE_ITEM = 'single item'
_SINGLE_ROW = 'single row'
_ROW_SLICE = 'row slice'
_DOUBLE_SLICE = 'double slice'

_SMALLEST_BLOCK = 1

_SLEEP_TIME = 0.2


class InvalidColorSpec(Exception):
    """
    Error for a color value that is not a number.

    """
    pass


def show_color(red, green, blue):
    """
    Show a given color in the IPython Notebook.

    Parameters
    ----------
    red, green, blue : int
        Integers on the range [0 - 255].

    """
    div = ('<div style="height: 60px; min-width: 200px; '
           'background-color: {0}"></div>')
    display(HTML(div.format(_RGB.format(red, green, blue))))


def embed_colorpicker():
    """
    Embed the web page www.colorpicker.com inside the IPython Notebook.

    """
    iframe = ('<iframe src="http://www.colorpicker.com/" '
              'width="100%" height="550px"></iframe>')
    display(HTML(iframe))


class Block(object):
    """
    A colored square.

    Parameters
    ----------
    red, green, blue : int
        Integers on the range [0 - 255].
    size : int, optional
        Length of the sides of this block in pixels. One is the lower limit.

    Attributes
    ----------
    red, green, blue : int
        The color values for this `Block`. The color of the `Block` can be
        updated by assigning new values to these attributes.
    rgb : tuple of int
        Tuple of (red, green, blue) values. Can be used to set all the colors
        at once.
    row, col : int
        The zero-based grid position of this `Block`.
    size : int
        Length of the sides of this block in pixels. The block size can be
        changed by modifying this attribute. Note that one is the lower limit.

    """

    def __init__(self, red, green, blue, size=20):
        self.red = red
        self.green = green
        self.blue = blue
        self.size = size

        self._row = None
        self._col = None

    @staticmethod
    def _check_value(value):
        """
        Check that a value is a number and constrain it to [0 - 255].

        """
        if not isinstance(value, numbers.Number):
            s = 'value must be a number. got {0}.'.format(value)
            raise InvalidColorSpec(s)

        return int(round(min(255, max(0, value))))

    @property
    def red(self):
        return self._red

    @red.setter
    def red(self, value):
        value = self._check_value(value)
        self._red = value

    @property
    def green(self):
        return self._green

    @green.setter
    def green(self, value):
        value = self._check_value(value)
        self._green = value

    @property
    def blue(self):
        return self._blue

    @blue.setter
    def blue(self, value):
        value = self._check_value(value)
        self._blue = value

    @property
    def rgb(self):
        return (self._red, self._green, self._blue)

    @rgb.setter
    def rgb(self, colors):
        if len(colors) != 3:
            s = 'Setting colors requires three values: (red, green, blue).'
            raise ValueError(s)

        self.red, self.green, self.blue = colors

    @property
    def row(self):
        return self._row

    @property
    def col(self):
        return self._col

    @property
    def size(self):
        return self._size

    @size.setter
    def size(self, size):
        self._size = max(_SMALLEST_BLOCK, size)

    def set_colors(self, red, green, blue):
        """
        Updated block colors.

        Parameters
        ----------
        red, green, blue : int
            Integers on the range [0 - 255].

        """
        self.red = red
        self.green = green
        self.blue = blue

    @property
    def _td(self):
        """
        The HTML for a table cell with the background color of this Block.

        """
        title = _TITLE.format(self._row, self._col,
                              self._red, self._green, self._blue)
        rgb = _RGB.format(self._red, self._green, self._blue)
        return _TD.format(title, self._size, rgb)

    def _repr_html_(self):
        return _TABLE.format(uuid.uuid4(), 0, _TR.format(self._td))

    def show(self):
        display(HTML(self._repr_html_()))

    def __str__(self):
        s = ['{0}'.format(self.__class__.__name__),
             'Color: ({0}, {1}, {2})'.format(self._red,
                                             self._green,
                                             self._blue)]

        # add position information if we have it
        if self._row is not None:
            s[0] += ' [{0}, {1}]'.format(self._row, self._col)

        return os.linesep.join(s)


class BlockGrid(object):
    """
    A grid of blocks whose colors can be individually controlled.

    Parameters
    ----------
    width : int
        Number of blocks wide to make the grid.
    height : int
        Number of blocks high to make the grid.
    fill : tuple of int, optional
        An optional initial color for the grid, defaults to black.
        Specified as a tuple of (red, green, blue). E.g.: (10, 234, 198)
    block_size : int, optional
        Length of the sides of grid blocks in pixels. One is the lower limit.
    lines_on : bool, optional
        Whether or not to display lines between blocks.

    Attributes
    ----------
    width : int
        Number of blocks along the width of the grid.
    height : int
        Number of blocks along the height of the grid.
    shape : tuple of int
        A tuple of (width, height).
    block_size : int
        Length of the sides of grid blocks in pixels. The block size can be
        changed by modifying this attribute. Note that one is the lower limit.
    lines_on : bool
        Whether lines are shown between blocks when the grid is displayed.
        This attribute can used to toggle the whether the lines appear.

    """

    def __init__(self, width, height, fill=(0, 0, 0),
                 block_size=20, lines_on=True):
        self._width = width
        self._height = height
        self._block_size = block_size
        self.lines_on = lines_on
        self._initialize_grid(fill)

    def _initialize_grid(self, fill):
        grid = [[Block(*fill, size=self._block_size)
                for col in xrange(self.width)]
                for row in xrange(self.height)]

        self._grid = grid

    @property
    def width(self):
        return self._width

    @property
    def height(self):
        return self._height

    @property
    def shape(self):
        return (self._width, self._height)

    @property
    def block_size(self):
        return self._block_size

    @block_size.setter
    def block_size(self, size):
        self._block_size = size

        for block in self:
            block.size = size

    @property
    def lines_on(self):
        return self._lines_on

    @lines_on.setter
    def lines_on(self, value):
        if value not in (0, 1):
            s = 'lines_on may only be True or False.'
            raise ValueError(s)

        self._lines_on = value

    @classmethod
    def _view_from_grid(cls, grid):
        """
        Make a new grid from a list of lists of Block objects.

        """
        new_width = len(grid[0])
        new_height = len(grid)

        new_BG = cls(new_width, new_height)
        new_BG._grid = grid

        return new_BG

    @staticmethod
    def _categorize_index(index):
        """
        Used by __getitem__ and __setitem__ to determine whether the user
        is asking for a single item, single row, or some kind of slice.

        """
        if isinstance(index, int):
            return _SINGLE_ROW

        elif isinstance(index, slice):
            return _ROW_SLICE

        elif isinstance(index, tuple):
            if len(index) > 2:
                s = 'Invalid index, too many dimensions.'
                raise IndexError(s)

            elif len(index) == 1:
                s = 'Single indices must be integers, not tuple.'
                raise TypeError(s)

            if isinstance(index[0], slice):
                if isinstance(index[1], (int, slice)):
                    return _DOUBLE_SLICE

            if isinstance(index[1], slice):
                if isinstance(index[0], (int, slice)):
                    return _DOUBLE_SLICE

            elif isinstance(index[0], int) and isinstance(index[0], int):
                return _SINGLE_ITEM

        raise IndexError('Invalid index.')

    def __getitem__(self, index):
        ind_cat = self._categorize_index(index)

        if ind_cat == _SINGLE_ROW:
            return self._view_from_grid([self._grid[index]])

        elif ind_cat == _SINGLE_ITEM:
            block = self._grid[index[0]][index[1]]
            block._row, block._col = index
            return block

        elif ind_cat == _ROW_SLICE:
            return self._view_from_grid(self._grid[index])

        elif ind_cat == _DOUBLE_SLICE:
            new_grid = self._get_double_slice(index)
            return self._view_from_grid(new_grid)

    def __setitem__(self, index, value):
        if len(value) != 3:
            s = 'Assigned value must have three integers. got {0}.'
            raise ValueError(s.format(value))

        ind_cat = self._categorize_index(index)

        if ind_cat == _SINGLE_ROW:
            map(lambda b: b.set_colors(*value), self._grid[index])

        elif ind_cat == _SINGLE_ITEM:
            self._grid[index[0]][index[1]].set_colors(*value)

        else:
            if ind_cat == _ROW_SLICE:
                sub_grid = self._grid[index]

            elif ind_cat == _DOUBLE_SLICE:
                sub_grid = self._get_double_slice(index)

            map(lambda b: b.set_colors(*value), itertools.chain(*sub_grid))

    def _get_double_slice(self, index):
        sl_height, sl_width = index

        if isinstance(sl_width, int):
            if sl_width == -1:
                sl_width = slice(sl_width, None)
            else:
                sl_width = slice(sl_width, sl_width + 1)

        if isinstance(sl_height, int):
            if sl_height == -1:
                sl_height = slice(sl_height, None)
            else:
                sl_height = slice(sl_height, sl_height + 1)

        rows = self._grid[sl_height]
        grid = [r[sl_width] for r in rows]

        return grid

    def __iter__(self):
        for r in xrange(self.height):
            for c in xrange(self.width):
                yield self[r, c]

    @property
    def animate(self):
        """
        Iterate over this property to have your changes to the grid
        animated in the IPython Notebook.

        """
        for block in self:
            self.show()
            time.sleep(_SLEEP_TIME)
            yield block
            clear_output()
        self.show()

    def _repr_html_(self):
        rows = range(self._height)
        cols = range(self._width)

        html = reduce(iadd,
                      (_TR.format(reduce(iadd,
                                         (self[r, c]._td
                                          for c in cols)))
                       for r in rows))

        return _TABLE.format(uuid.uuid4(), int(self._lines_on), html)

    def __str__(self):
        s = ['{0}'.format(self.__class__.__name__),
             'Shape: {0}'.format(self.shape)]

        return os.linesep.join(s)

    def copy(self):
        """
        Returns an independent copy of this BlockGrid.

        """
        return copy.deepcopy(self)

    def show(self):
        """
        Display colored grid as an HTML table.

        """
        display(HTML(self._repr_html_()))

    def flash(self):
        """
        Display the grid for a short time. Useful for making an animation.

        """
        self.show()
        time.sleep(_SLEEP_TIME)
        clear_output()

    def to_text(self, filename=None):
        """
        Write a text file containing the size and block color information
        for this grid.

        If no file name is given the text is sent to stdout.

        Parameters
        ----------
        filename : str, optional
            File into which data will be written. Will be overwritten if
            it already exists.

        """
        if filename:
            f = open(filename, 'w')
        else:
            f = sys.stdout

        s = ['# width height', '{0} {1}'.format(self.width, self.height),
             '# block size', '{0}'.format(self.block_size),
             '# initial color', '0 0 0',
             '# row column red green blue']
        f.write(os.linesep.join(s) + os.linesep)

        for block in self:
            things = [str(x) for x in (block.row, block.col) + block.rgb]
            f.write(' '.join(things) + os.linesep)

        if filename:
            f.close()


class Pixel(Block):
    @property
    def x(self):
        """
        Horizontal coordinate of Pixel.

        """
        return self._col

    @property
    def y(self):
        """
        Vertical coordinate of Pixel.

        """
        return self._row

    @property
    def _td(self):
        """
        The HTML for a table cell with the background color of this Pixel.

        """
        title = _TITLE.format(self._col, self._row,
                              self._red, self._green, self._blue)
        rgb = _RGB.format(self._red, self._green, self._blue)
        return _TD.format(title, self._size, rgb)

    def __str__(self):
        s = ['{0}'.format(self.__class__.__name__),
             'Color: ({0}, {1}, {2})'.format(self._red,
                                             self._green,
                                             self._blue)]

        # add position information if we have it
        if self._row is not None:
            s[0] += ' [{0}, {1}]'.format(self._col, self._row)

        return os.linesep.join(s)


class ImageGrid(BlockGrid):
    """
    A grid of blocks whose colors can be individually controlled.

    Parameters
    ----------
    width : int
        Number of blocks wide to make the grid.
    height : int
        Number of blocks high to make the grid.
    fill : tuple of int, optional
        An optional initial color for the grid, defaults to black.
        Specified as a tuple of (red, green, blue). E.g.: (10, 234, 198)
    block_size : int, optional
        Length of the sides of grid blocks in pixels. One is the lower limit.
    lines_on : bool, optional
        Whether or not to display lines between blocks.
    origin : {'lower-left', 'upper-left'}
        Set the location of the grid origin.

    Attributes
    ----------
    width : int
        Number of blocks along the width of the grid.
    height : int
        Number of blocks along the height of the grid.
    shape : tuple of int
        A tuple of (width, height).
    block_size : int
        Length of the sides of grid blocks in pixels.
    lines_on : bool
        Whether lines are shown between blocks when the grid is displayed.
        This attribute can used to toggle the whether the lines appear.
    origin : str
        The location of the grid origin.

    """

    def __init__(self, width, height, fill=(0, 0, 0),
                 block_size=20, lines_on=True, origin='lower-left'):
        super(ImageGrid, self).__init__(width, height, fill,
                                        block_size, lines_on)

        if origin not in ('lower-left', 'upper-left'):
            s = "origin keyword must be one of {'lower-left', 'upper-left'}."
            raise ValueError(s)

        self._origin = origin

    def _initialize_grid(self, fill):
        grid = [[Pixel(*fill, size=self._block_size)
                for col in xrange(self.width)]
                for row in xrange(self.height)]

        self._grid = grid

    @property
    def block_size(self):
        return self._block_size

    @property
    def origin(self):
        return self._origin

    def _transform_index(self, index):
        """
        Transform a single-item index from Python style coordinates to
        image style coordinates in which the first item refers to column and
        the second item refers to row. Also takes into account the
        location of the origin.

        """
        # in ImageGrid index is guaranteed to be a tuple.

        # first thing, switch the coordinates since ImageGrid is column
        # major and ._grid is row major.
        new_ind = [index[1], index[0]]

        # now take into account that the ImageGrid origin may be lower-left,
        # while the ._grid origin is upper-left.
        if self._origin == 'lower-left':
            new_ind[0] = self._height - new_ind[0] - 1

        return tuple(new_ind)

    def __getitem__(self, index):
        ind_cat = self._categorize_index(index)

        # ImageGrid will only support single item indexing and 2D slices
        if ind_cat not in (_DOUBLE_SLICE, _SINGLE_ITEM):
            s = 'ImageGrid only supports 2D indexing.'
            raise IndexError(s)

        if ind_cat == _SINGLE_ITEM:
            real_index = self._transform_index(index)
            pixel = self._grid[real_index[0]][real_index[1]]
            pixel._col, pixel._row = index
            return pixel

        elif ind_cat == _DOUBLE_SLICE:
            new_grid = self._get_double_slice(index)
            return self._view_from_grid(new_grid)

    def __setitem__(self, index, value):
        if len(value) != 3:
            s = 'Assigned value must have three integers. got {0}.'
            raise ValueError(s.format(value))

        pixels = self[index]

        if isinstance(pixels, Pixel):
            pixels.set_colors(*value)

        else:
            map(lambda p: p.set_colors(*value), itertools.chain(*pixels._grid))

    def _get_double_slice(self, index):
        cslice, rslice = index

        if isinstance(rslice, int):
            if rslice == -1:
                rslice = slice(rslice, None)
            else:
                rslice = slice(rslice, rslice + 1)

        if isinstance(cslice, int):
            if cslice == -1:
                cslice = slice(cslice, None)
            else:
                cslice = slice(cslice, cslice + 1)

        rows = range(self._height)[rslice]
        if self._origin == 'lower-left':
            rows = rows[::-1]

        cols = range(self._width)[cslice]

        new_grid = [[self[c, r] for c in cols] for r in rows]

        return new_grid

    def __iter__(self):
        for col in xrange(self.width):
            for row in xrange(self.height):
                yield self[col, row]

    def _repr_html_(self):
        rows = range(self._height)
        cols = range(self._width)

        if self._origin == 'lower-left':
            rows = rows[::-1]

        html = reduce(iadd,
                      (_TR.format(reduce(iadd,
                                         (self[c, r]._td
                                          for c in cols)))
                       for r in rows))

        return _TABLE.format(uuid.uuid4(), int(self._lines_on), html)


# As a convenience, provide the named HTML colors as a dictionary.
colors = \
    {'AliceBlue': (240, 248, 255),
     'AntiqueWhite': (250, 235, 215),
     'Aqua': (0, 255, 255),
     'Aquamarine': (127, 255, 212),
     'Azure': (240, 255, 255),
     'Beige': (245, 245, 220),
     'Bisque': (255, 228, 196),
     'Black': (0, 0, 0),
     'BlanchedAlmond': (255, 235, 205),
     'Blue': (0, 0, 255),
     'BlueViolet': (138, 43, 226),
     'Brown': (165, 42, 42),
     'BurlyWood': (222, 184, 135),
     'CadetBlue': (95, 158, 160),
     'Chartreuse': (127, 255, 0),
     'Chocolate': (210, 105, 30),
     'Coral': (255, 127, 80),
     'CornflowerBlue': (100, 149, 237),
     'Cornsilk': (255, 248, 220),
     'Crimson': (220, 20, 60),
     'Cyan': (0, 255, 255),
     'DarkBlue': (0, 0, 139),
     'DarkCyan': (0, 139, 139),
     'DarkGoldenrod': (184, 134, 11),
     'DarkGray': (169, 169, 169),
     'DarkGreen': (0, 100, 0),
     'DarkKhaki': (189, 183, 107),
     'DarkMagenta': (139, 0, 139),
     'DarkOliveGreen': (85, 107, 47),
     'DarkOrange': (255, 140, 0),
     'DarkOrchid': (153, 50, 204),
     'DarkRed': (139, 0, 0),
     'DarkSalmon': (233, 150, 122),
     'DarkSeaGreen': (143, 188, 143),
     'DarkSlateBlue': (72, 61, 139),
     'DarkSlateGray': (47, 79, 79),
     'DarkTurquoise': (0, 206, 209),
     'DarkViolet': (148, 0, 211),
     'DeepPink': (255, 20, 147),
     'DeepSkyBlue': (0, 191, 255),
     'DimGray': (105, 105, 105),
     'DodgerBlue': (30, 144, 255),
     'FireBrick': (178, 34, 34),
     'FloralWhite': (255, 250, 240),
     'ForestGreen': (34, 139, 34),
     'Fuchsia': (255, 0, 255),
     'Gainsboro': (220, 220, 220),
     'GhostWhite': (248, 248, 255),
     'Gold': (255, 215, 0),
     'Goldenrod': (218, 165, 32),
     'Gray': (128, 128, 128),
     'Green': (0, 128, 0),
     'GreenYellow': (173, 255, 47),
     'Honeydew': (240, 255, 240),
     'HotPink': (255, 105, 180),
     'IndianRed': (205, 92, 92),
     'Indigo': (75, 0, 130),
     'Ivory': (255, 255, 240),
     'Khaki': (240, 230, 140),
     'Lavender': (230, 230, 250),
     'LavenderBlush': (255, 240, 245),
     'LawnGreen': (124, 252, 0),
     'LemonChiffon': (255, 250, 205),
     'LightBlue': (173, 216, 230),
     'LightCoral': (240, 128, 128),
     'LightCyan': (224, 255, 255),
     'LightGoldenrodYellow': (250, 250, 210),
     'LightGray': (211, 211, 211),
     'LightGreen': (144, 238, 144),
     'LightPink': (255, 182, 193),
     'LightSalmon': (255, 160, 122),
     'LightSeaGreen': (32, 178, 170),
     'LightSkyBlue': (135, 206, 250),
     'LightSlateGray': (119, 136, 153),
     'LightSteelBlue': (176, 196, 222),
     'LightYellow': (255, 255, 224),
     'Lime': (0, 255, 0),
     'LimeGreen': (50, 205, 50),
     'Linen': (250, 240, 230),
     'Magenta': (255, 0, 255),
     'Maroon': (128, 0, 0),
     'MediumAquamarine': (102, 205, 170),
     'MediumBlue': (0, 0, 205),
     'MediumOrchid': (186, 85, 211),
     'MediumPurple': (147, 112, 219),
     'MediumSeaGreen': (60, 179, 113),
     'MediumSlateBlue': (123, 104, 238),
     'MediumSpringGreen': (0, 250, 154),
     'MediumTurquoise': (72, 209, 204),
     'MediumVioletRed': (199, 21, 133),
     'MidnightBlue': (25, 25, 112),
     'MintCream': (245, 255, 250),
     'MistyRose': (255, 228, 225),
     'Moccasin': (255, 228, 181),
     'NavajoWhite': (255, 222, 173),
     'Navy': (0, 0, 128),
     'OldLace': (253, 245, 230),
     'Olive': (128, 128, 0),
     'OliveDrab': (107, 142, 35),
     'Orange': (255, 165, 0),
     'OrangeRed': (255, 69, 0),
     'Orchid': (218, 112, 214),
     'PaleGoldenrod': (238, 232, 170),
     'PaleGreen': (152, 251, 152),
     'PaleTurquoise': (175, 238, 238),
     'PaleVioletRed': (219, 112, 147),
     'PapayaWhip': (255, 239, 213),
     'PeachPuff': (255, 218, 185),
     'Peru': (205, 133, 63),
     'Pink': (255, 192, 203),
     'Plum': (221, 160, 221),
     'PowderBlue': (176, 224, 230),
     'Purple': (128, 0, 128),
     'Red': (255, 0, 0),
     'RosyBrown': (188, 143, 143),
     'RoyalBlue': (65, 105, 225),
     'SaddleBrown': (139, 69, 19),
     'Salmon': (250, 128, 114),
     'SandyBrown': (244, 164, 96),
     'SeaGreen': (46, 139, 87),
     'Seashell': (255, 245, 238),
     'Sienna': (160, 82, 45),
     'Silver': (192, 192, 192),
     'SkyBlue': (135, 206, 235),
     'SlateBlue': (106, 90, 205),
     'SlateGray': (112, 128, 144),
     'Snow': (255, 250, 250),
     'SpringGreen': (0, 255, 127),
     'SteelBlue': (70, 130, 180),
     'Tan': (210, 180, 140),
     'Teal': (0, 128, 128),
     'Thistle': (216, 191, 216),
     'Tomato': (255, 99, 71),
     'Turquoise': (64, 224, 208),
     'Violet': (238, 130, 238),
     'Wheat': (245, 222, 179),
     'White': (255, 255, 255),
     'WhiteSmoke': (245, 245, 245),
     'Yellow': (255, 255, 0),
     'YellowGreen': (154, 205, 50)}
