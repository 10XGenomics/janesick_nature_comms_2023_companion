import math
from collections import namedtuple, Counter

import numpy as np
import scanpy as sc
import pandas as pd

## Image manipulation and geometry
from cv2 import perspectiveTransform
from ome_types import from_tiff
from skimage.measure import points_in_poly

from scipy.spatial import KDTree, ConvexHull
from scipy.sparse import csc_matrix
from shapely.geometry import LinearRing

## Plotting imports
from matplotlib import pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.colors import LinearSegmentedColormap, to_hex, Colormap
from matplotlib.axes import Axes


def hexlist_to_cmap(name: str, hex_list: list[str]) -> Colormap:
    """
    Given a list of hex colors produces a matplotlib discrete colormap.

    """
    colors = [to_hex(hex_color) for hex_color in hex_list]
    cmap = LinearSegmentedColormap.from_list(name=name, colors=colors, N=len(colors))
    return cmap


def polygons_coords_to_patch_collection(
    poly_coords: list[list[float]], **patch_collection_kwargs
) -> PatchCollection:
    """
    Generates a matplotlib PatchCollection from a set of polygons. A helper function for plot_polygons.
    """
    patches = [Polygon(hexagon) for hexagon in poly_coords]
    collection = PatchCollection(patches, **patch_collection_kwargs)
    return collection


def plot_polygons(
    poly_coords: list[list[float]],
    ax: Axes = None,
    auto_set_lim: bool = True,
    **patch_collection_kwargs
):
    """
    Plots a set of polygons. Polygons may be styled using patch_collection_kwargs passed to Collection.
    e.g.
    array, facecolor, edgecolor, linewidth....
    https://matplotlib.org/stable/api/collections_api.html#matplotlib.collections.Collection

    pitfall:
    when auto_set_lim == False, a plot will no rescale to include polygons by default.
    """
    if ax is None:
        ax = plt.gca()

    ax.add_collection(
        polygons_coords_to_patch_collection(poly_coords, **patch_collection_kwargs)
    )
    if auto_set_lim:
        ax.autoscale_view()


def hex_corner_offset(corner: int, size: float) -> np.ndarray[float, (2,)]:
    """
    Helper function for polygon_corners, calculates the offset of each hexagon vertex from the center.
    """
    angle_deg = 60 * corner - 30
    angle_rad = math.pi / 180 * angle_deg
    offset = np.array([size[0] * math.cos(angle_rad), size[1] * math.sin(angle_rad)])
    return offset


def polygon_corners(
    center: np.ndarray[float, (2,)], size: float
) -> np.ndarray[float, (6, 2)]:
    """
    Generates the polygon vertices for a hexegon of size at origin center.
    Used for plotting Visium space filling Hexagons.
    Reference for working with Hexagons: https://www.redblobgames.com/grids/hexagons/
    Visium spots are organized such that space filling hexagons are POINTY.
    """
    return np.array([center + hex_corner_offset(i, size) for i in range(0, 6)])


## The hexcodes for the colors used within the paper.
## The celltypes are mapped to the colors for plotting convenience.
celltypes, hex_codes = np.array(
    sorted(
        [
            [
                "DCIS_1",
                "#fb34cd",
            ],
            [
                "DCIS_2",
                "#FE664D",
            ],
            [
                "Myoepi_ACTA2+",
                "#009203",
            ],
            [
                "Myoepi_KRT15+",
                "#66c102",
            ],
            [
                "Invasive_Tumor",
                "#ff002a",
            ],
            [
                "Prolif_Invasive_Tumor",
                "#8e0119",
            ],
            [
                "T_Cell_&_Tumor_Hybrid",
                "#CB6035",
            ],
            [
                "Stromal",
                "#e5e022",
            ],
            [
                "Stromal_&_T_Cell_Hybrid",
                "#C1A029",
            ],
            [
                "CD4+_T_Cells",
                "#4fa9ff",
            ],
            [
                "CD8+_T_Cells",
                "#1068be",
            ],
            [
                "B_Cells",
                "#565DFD",
            ],
            [
                "Macrophages_1",
                "#10686f",
            ],
            [
                "Macrophages_2",
                "#3694a8",
            ],
            [
                "IRF7+_DCs",
                "#9f50f9",
            ],
            [
                "LAMP3+_DCs",
                "#AB76AE",
            ],
            [
                "Mast_Cells",
                "#999999",
            ],
            [
                "Perivascular-Like",
                "#515151",
            ],
            ["Endothelial", "#01257b"],
            ["Unlabeled", "#ffa5aa"],
        ],
        key=lambda x: x[0],
    )
).T
## used for converting list of celltypes to categorical.
ctype_to_code_map = {ctype: i for i, ctype in enumerate(celltypes)}
## used for converting list of celltypes to colors.
ctype_hex_map = dict(zip(celltypes, hex_codes))
## matplotlib colormap
ctype_cmap = hexlist_to_cmap("celltype_cmap", hex_codes)


def unique_encode(
    array: np.ndarray[str], return_inverse: bool = True
) -> (np.ndarray[str], np.ndarray[int]):
    """np.unique is slower then list(set(array)) thus this is an emulator function of np.unique.
    for fast encoding of categorical arrays.
    """
    unique_values = np.array(list(set(array)))
    if not return_inverse:
        return unique_values
    inverse_map = {val: i for i, val in enumerate(unique_values)}
    inverse_codes = np.array([inverse_map[i] for i in array])
    return (unique_values, inverse_codes)


def get_xenium_to_morphology_transform_from_xenium_morphology_image(
    morphology_path: str, pyramid_level: int = 0
) -> np.ndarray[np.float32]:
    """
    The Xenium coordinate system has its origin at (0,0) set to the top left of the Xenium morphology image.
    The image coordinate system is in uM and must be converted to px. This produces a transformation matrix for
    that conversion.

    https://kb.10xgenomics.com/hc/en-us/articles/11636252598925-What-are-the-Xenium-image-scale-factors-
    """
    morphology_metadata = from_tiff(morphology_path)
    origin_x = morphology_metadata.plates[0].well_origin_x
    origin_y = morphology_metadata.plates[0].well_origin_y
    physical_size_x = morphology_metadata.images[0].pixels.physical_size_x
    physical_size_y = morphology_metadata.images[0].pixels.physical_size_y
    ##Homography matrix describing scaling by pixel-size and translating by offset.
    return np.array(
        [
            [1 / physical_size_x / 2**pyramid_level, 0, 0],
            [0, 1 / physical_size_y / 2**pyramid_level, 0],
            [0, 0, 1],
        ]
    )


def get_xenium_capture_polygon_um(
    path_to_xenium_bundle: str,
) -> np.ndarray[int, (4, 2)]:
    """Produces a polygon which describes the Xenium morphology image shape.
    Polygons are easily transformable. Thus, this polygon is useful for sub-setting coordinates
    that lie within an image under arbitrary transformation, and of course for debug purposes.
    """
    morphology_metadata = from_tiff(path_to_xenium_bundle / "morphology_mip.ome.tif")
    image_shape_px = np.array(
        [
            morphology_metadata.images[0].pixels.size_x,
            morphology_metadata.images[0].pixels.size_y,
        ]
    )
    image_shape_um = (
        image_shape_px * morphology_metadata.images[0].pixels.physical_size_x
    )
    xenium_capture_polygon_um = np.array(
        [(0, 0), (0, image_shape_um[1]), image_shape_um, (image_shape_um[0], 0)]
    )
    return xenium_capture_polygon_um


def transform_coordinates(
    coords: np.ndarray[np.float32, (None, 2)], transform_matrix: np.ndarray[np.float32]
) -> np.ndarray[np.float32, (None, 2)]:
    """Transforms coordinates using transform.
    https://en.wikipedia.org/wiki/Homography_(computer_vision)
    """
    return perspectiveTransform(
        coords.reshape(-1, 1, 2).astype(np.float32), transform_matrix
    )[:, 0, :]


def get_median_spot_to_spot_distance_from_centroids(
    centroids: np.ndarray[np.float32, (None, 2)]
) -> float:
    """
    Calculates the median distance between all coordinates and their closest neighbor.
    Usefull for quickly determining the median spot <-> spot distance for Visium spots in the current coordinate system.
    Can also be directly calculated using the full resolution image pixel size and the theoretical spot <-> spot distance.
    This distance is used to generate a rough estimate of the Visium capture area and to produce space filling
    visium spot Hexagons.
    """
    tree = KDTree(centroids)
    distances, neighbor = tree.query(centroids, k=2)
    return np.median(distances[0, 1])


def generate_space_filling_visium_polygons(
    spot_coords: np.ndarray[np.float32, (None, 2)], spot_to_spot_distance: float = None
) -> list[np.ndarray[float, (6, 2)]]:
    """
    Generates hexagon polygons for every visium spot, by default the polygons are fully space filling.
    That is to say the size of the hexagon is set using the spot <-> spot distance.
    These polygons are for display purposes only and should not be used for direct analysis.
    """
    if spot_to_spot_distance is None:
        spot_to_spot_distance = get_median_spot_to_spot_distance_from_centroids(
            spot_coords
        )
    size = spot_to_spot_distance / math.sqrt(3)
    visium_hexegon_polygons = np.array(
        [polygon_corners(coord, [size, size]) for coord in spot_coords]
    )
    return visium_hexegon_polygons


def get_visium_capture_polygon(
    spot_centroids: np.ndarray[np.float32, (None, 2)], spot_diameter: float
) -> np.ndarray[np.float32, (None, 2)]:
    """
    Generates a polygon representing the visium capture area; It calculates the convex hull of all Visium spots
    then pads the convex hall by the diameter of one spot. This is used for filtering data outside the capture area
    as well as visualization.

    https://en.wikipedia.org/wiki/Convex_hull
    """
    convex_hull = ConvexHull(spot_centroids)
    convex_poly = convex_hull.points[convex_hull.vertices]

    padded_convex_polly = np.array(
        LinearRing(convex_poly).buffer(spot_diameter).exterior.coords.xy
    ).T
    return padded_convex_polly


__OUTSIDE_VISIUM_CAPTURE_AREA_BARCODE__ = "Outside_Visium_Capture_Area"


def bin_xenium_data_to_visium_spots(
    visium_spot_centroids: np.ndarray[np.float32, (None, 2)],
    visium_spot_diameter: float,
    visium_barcodes: np.ndarray[str, (None, 2)],
    xenium_coords: np.ndarray[np.float32, (None, 2)],
    xenium_to_visium_transform: np.ndarray[np.float32],
) -> np.ndarray[str]:
    """
    Assigns xenium coordinates to visium spots using a transform into the visium coordinate system.
    Coordinates are assigned by nearest neighbor.
    e.g. each Xenium transcript/cell etc is assigned to the closest Visium barcode after transformation.
    Also filters coordinates that lie outside an estimated capture area.
    """
    xenium_coords_in_visium = transform_coordinates(
        xenium_coords, xenium_to_visium_transform
    )
    kdtree = KDTree(visium_spot_centroids)
    _, neighbors = kdtree.query(xenium_coords_in_visium, k=1)
    barcode_assignments = visium_barcodes[neighbors]

    visium_capture_polygon_full_res = get_visium_capture_polygon(
        spot_centroids=visium_spot_centroids, spot_diameter=visium_spot_diameter
    )
    barcode_assignments[
        ~points_in_poly(xenium_coords_in_visium, visium_capture_polygon_full_res)
    ] = __OUTSIDE_VISIUM_CAPTURE_AREA_BARCODE__
    return barcode_assignments


def generate_anndata_matrix_from_transcript_assignments(
    barcodes: np.ndarray[str],
    feature_names: np.ndarray[str],
) -> sc.AnnData:
    """Generates a count matrix as an anndata object from a set of barcodes and feature names.
    e.g. given transcript cell assignments and gene names it will produce a barcode feature matrix.
    """
    valid_feature_names, feature_codes = unique_encode(
        feature_names, return_inverse=True
    )
    valid_barcodes, barcode_codes = unique_encode(barcodes, return_inverse=True)

    mat_dict = Counter(zip(barcode_codes, feature_codes))
    row_ind_col_ind = np.array(list(mat_dict.keys()), dtype=int).T
    data = np.array(list(mat_dict.values()), dtype=np.float32)

    obs = pd.DataFrame(index=valid_barcodes)
    var = pd.DataFrame(index=valid_feature_names)
    adata = sc.AnnData(X=csc_matrix((data, row_ind_col_ind)), obs=obs, var=var)
    return adata
