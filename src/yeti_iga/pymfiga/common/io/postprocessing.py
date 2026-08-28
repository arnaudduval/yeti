from typing import Optional, Literal
import numpy as np
import logging
import os

logger = logging.getLogger("SRC")

try:
    from PIL import Image
except:
    logger.debug(
        "PIL module is not installed. Some functions in post-processing may be disabled"
    )


try:
    import pyvista as pv
except:
    logger.debug(
        "pyvista module is not installed. Some functions in post-processing may be disabled"
    )


class Postprocessing:
    @staticmethod
    def verify_folder_to_save(folder: Optional[str]):
        if folder is None:
            folder = os.path.join(os.getcwd(), "results/")
        if not os.path.isdir(folder):
            logger.info(f"Folder {folder} will be created")
            os.mkdir(folder)
        logger.info(f"File will be save at {folder}")
        return folder

    @staticmethod
    def crop_image(filename: str):
        im = np.array(Image.open(filename).convert("RGB"))
        colorY, colorX = np.where(np.all(im != [255, 255, 255], axis=2))
        top, bottom = min(colorY), max(colorY)
        left, right = min(colorX), max(colorX)
        Image.fromarray(im[top:bottom, left:right]).save(filename)

    @staticmethod
    def vtk2png(
        filename: str,
        fieldname: str,
        folder: Optional[str] = None,
        cmap: str = "viridis",
        format: Literal["vts", "vtk"] = "vtk",
        **kwargs,
    ):
        folder = Postprocessing.verify_folder_to_save(folder)
        filename_to_read = f"{folder}/{filename}.{format}"
        mesh = pv.read(filename_to_read)
        filename_to_save = f"{folder}/{filename}.png"
        mesh.plot(
            scalars=fieldname,
            cmap=cmap,
            screenshot=filename_to_save,
            off_screen=True,
            **kwargs,
        )
        Postprocessing.crop_image(filename_to_save)
