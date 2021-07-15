"""
Module notebook
================

A module with helper functions for Jupyterlab notebooks.

"""
import IPython
from IPython import display
from IPython.display import HTML
from IPython.display import Image


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def _src_from_data(data):
    """Base64 encodes image bytes for inclusion in an HTML img element.

    Very strongly inspired after:
    https://mindtrove.info/

    :param data: the image data
    :return: the src image code to embed in the HTML
    """
    img_obj = Image(data=data)
    for bundle in img_obj._repr_mimebundle_():
        for mimetype, b64value in bundle.items():
            if mimetype.startswith('image/'):
                return f'data:{mimetype};base64,{b64value}'


def gallery(images: list, row_height: str = 'auto') -> str:
    """Shows a set of images in a gallery that flexes with the width of the notebook.

    This works only with PNG images.

    Very strongly inspired after:
    https://mindtrove.info/

    :param images: list of URLs, bytes or IPython.Display.Image objects to display.
    :param row_height: CSS height value to assign to all images. Set to 'auto' by default to show images with their native dimensions. Set to a value like '250px' to make all rows in the gallery equal height.
    :return: the HTML code for displaying the gallery of images in Jupyterlab.
    """
    figures = []
    for image in images:
        if isinstance(image, IPython.core.display.Image):
            image = image.data
        if isinstance(image, bytes):
            src = _src_from_data(image)
            caption = ''
        else:
            src = image
            caption = f'<figcaption style="font-size: 0.6em">{image}</figcaption>'
        figures.append(f'''
            <figure style="margin: 5px !important;">
              <img src="{src}" style="height: {row_height}">
              {caption}
            </figure>
        ''')
    return HTML(data=f'''
        <div style="display: flex; flex-flow: row wrap; text-align: center;">
        {''.join(figures)}
        </div>
    ''')


def display_image_table(imgs, max_img_per_row=4):
    """Display a table filled with images.

    :param imgs: an iterable containing images
    :param max_img_per_row: the maximum number of images per row
    """

    table = """
<table>
    <tr>
"""
    for i, img in enumerate(imgs):
        if i % max_img_per_row == 0:
            table += """
            </tr>
            <tr>
            """
        table += f"<td>{img}</td>\n"
    table += """
    </tr>
</table>
"""
    return HTML(table)
