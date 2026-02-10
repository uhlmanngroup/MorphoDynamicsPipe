#uses view3 environment (no numba and napari=0.4.18)
import napari, os, skimage, pathlib
import numpy as np
import pandas as pd
from qtpy import QtWidgets
from qtpy.QtWidgets import QFileDialog, QWidget, QHBoxLayout, QVBoxLayout, QPushButton, QLineEdit, QLabel, QComboBox
from skimage.color import label2rgb
from magicgui import magicgui
from napari.types import LabelsData, ImageData, LayerDataTuple
from napari.layers import Layer, Labels, Image
from typing import List
import natsort
import dask_image.imread
import warnings
warnings.filterwarnings("ignore")

@magicgui(
    folder={"label": "Folder to open:\n1data/subfolder/", "mode": "d"},
    open_1_data={"widget_type": "CheckBox", "label": "1_data"},
    open_2_segmentation={"widget_type": "CheckBox", "label": "2_segmentation"},
    open_3b_tracking_images = {"widget_type": "CheckBox", "label": "3b_tracking_images"},
    open_3c_tracking_images_filtered = {"widget_type": "CheckBox", "label": "3c_tracking_images_filtered"},
    open_5_tracking_images_outlines = {"widget_type": "CheckBox", "label": "5_tracking_images_outlines"},
    lazy_loading={"widget_type": "CheckBox", "label": "Lazy loading"},
)
def open_folder(folder: pathlib.Path = '/nfs/research/uhlmann/', 
#                folder: pathlib.Path = '/Users/bwoodhams/Documents/plast_cell/1_data/',
                open_1_data:bool = True,
                open_2_segmentation:bool = False,
                open_3b_tracking_images:bool = False,
                open_3c_tracking_images_filtered:bool = False,
                open_5_tracking_images_outlines:bool = False,
                lazy_loading:bool = True,
                ) -> LayerDataTuple: #if put List here it doesn't work
    list_of_files = [each for each in os.listdir(folder) if each.endswith('.tif')]
    list_of_files = natsort.natsorted(list_of_files)

    dict_image_or_labels = {'1_data':'image', 
                            '2_segmentation':'labels', 
                            '3b_tracking_images':'labels', 
                            '3c_tracking_images_filtered':'labels', 
                            '5_tracking_images_outlines':'labels'}
    output_list = []

    if lazy_loading:
        for open_this, name in zip([open_1_data, open_2_segmentation, open_3b_tracking_images, open_3c_tracking_images_filtered, open_5_tracking_images_outlines], 
                                ['1_data', '2_segmentation', '3b_tracking_images', '3c_tracking_images_filtered', '5_tracking_images_outlines']):
            if open_this:
                this_folder = str(folder).replace('1_data', name)
                mylistofimages = dask_image.imread.imread(os.path.join(this_folder, '*.tif'))
                output_list.append((mylistofimages, {'name': name}, dict_image_or_labels[name]))
    else:
        for open_this, name in zip([open_1_data, open_2_segmentation, open_3b_tracking_images, open_3c_tracking_images_filtered, open_5_tracking_images_outlines], 
                                ['1_data', '2_segmentation', '3b_tracking_images', '3c_tracking_images_filtered', '5_tracking_images_outlines']):
            if open_this:
                this_folder = str(folder).replace('1_data', name)
                mylistofimages = np.array([skimage.io.imread(os.path.join(this_folder,file)) for file in list_of_files])
                output_list.append((mylistofimages, {'name': name}, dict_image_or_labels[name]))

    return output_list

viewer = napari.Viewer()
viewer.window.add_dock_widget(open_folder, area='bottom')

# --- new: CSV-to-layer dock widget -----------------------------------------
def _create_csv_to_layer_widget(viewer):
    w = QWidget()
    layout = QVBoxLayout()
    # file chooser row
    row = QHBoxLayout()
    row.addWidget(QLabel("CSV:"))
    file_edit = QLineEdit()
    file_edit.setReadOnly(True)
    row.addWidget(file_edit)
    browse_btn = QPushButton("Browse")
    row.addWidget(browse_btn)
    layout.addLayout(row)
    # columns row
    col_row = QHBoxLayout()
    col_row.addWidget(QLabel("Column:"))
    column_combo = QComboBox()
    column_combo.setEditable(False)
    col_row.addWidget(column_combo)
    layout.addLayout(col_row)
    # action
    add_btn = QPushButton("Add CSV column as image layer")
    layout.addWidget(add_btn)
    w.setLayout(layout)

    # internal storage
    w._df = None
    w._csv_path = None

    def on_browse():
        path, _ = QFileDialog.getOpenFileName(w, "Select CSV", "", "CSV Files (*.csv);;All Files (*)")
        if not path:
            return
        try:
            df = pd.read_csv(path)
        except Exception as e:
            print("Failed to read CSV:", e)
            return
        w._df = df
        w._csv_path = path
        file_edit.setText(path)
        column_combo.clear()
        column_combo.addItems(list(df.columns.astype(str)))

    def on_add():
        if w._df is None:
            print("No CSV loaded")
            return
        col = column_combo.currentText()
        if not col:
            print("No column selected")
            return
        # find labels layer named '3c_tracking_images_filtered'
        if '3c_tracking_images_filtered' not in viewer.layers:
            print("Layer '3c_tracking_images_filtered' not found")
            return
        labels_layer = viewer.layers['3c_tracking_images_filtered']
        labels = np.array(labels_layer.data)  # bring into memory
        df = w._df

        # determine id and frame columns
        id_col = 'cellID' if 'cellID' in df.columns else df.columns[0]
        frame_col = 'frame_id_T' if 'frame_id_T' in df.columns else None

        # time-aware mapping
        try:
            if frame_col is None:
                # no time column: behave as single-frame mapping
                keys = df[id_col].astype(int).values
                values = df[col].values
                mapping = dict(zip(keys, values))
                out = np.zeros_like(labels, dtype=float)
                for lab in np.unique(labels):
                    if lab == 0:
                        continue
                    val = mapping.get(int(lab), 0.0)
                    if val != 0:
                        out[labels == lab] = float(val)
            else:
                # ensure labels has a time axis; if 2D add a time axis of length 1
                if labels.ndim == 2:
                    labels_t = labels[np.newaxis, ...]
                    single_frame = True
                else:
                    labels_t = labels
                    single_frame = False
                n_frames = labels_t.shape[0]

                # prepare dataframe subset and cast to ints
                df_sub = df[[id_col, frame_col, col]].dropna(subset=[id_col, frame_col])
                df_sub[id_col] = df_sub[id_col].astype(int)
                df_sub[frame_col] = df_sub[frame_col].astype(int)

                # adjust 1-based frame indices if needed
                min_f = int(df_sub[frame_col].min())
                max_f = int(df_sub[frame_col].max())
                if max_f >= n_frames and min_f > 0:
                    df_sub[frame_col] = df_sub[frame_col] - 1

                # keep only rows that map into available frames
                df_sub = df_sub[(df_sub[frame_col] >= 0) & (df_sub[frame_col] < n_frames)]

                out = np.zeros_like(labels_t, dtype=float)
                # iterate CSV rows and fill per-frame masks
                for _, row in df_sub.iterrows():
                    f = int(row[frame_col])
                    labid = int(row[id_col])
                    try:
                        val = float(row[col])
                    except Exception:
                        val = row[col]
                    if val != 0:
                        out[f][labels_t[f] == labid] = float(val)

                if single_frame:
                    out = out[0]

            viewer.add_image(out, name=f"{col}_from_csv", blending='additive', colormap='magma')
        except Exception as e:
            print("Failed to create time-aware layer:", e)

    browse_btn.clicked.connect(on_browse)
    add_btn.clicked.connect(on_add)
    return w

csv_widget = _create_csv_to_layer_widget(viewer)
viewer.window.add_dock_widget(csv_widget, area='right')
# ---------------------------------------------------------------------------

napari.run()
