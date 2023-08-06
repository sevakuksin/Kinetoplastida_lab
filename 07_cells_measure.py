# import os
# import sys

import matplotlib.pyplot as plt
# import matplotlib.patches as patches
import matplotlib.patches as mpatches
# from matplotlib import colors
# from matplotlib.colorsz import LinearSegmentedColormap

import numpy as np
from numpy import ma

from skimage import io
from skimage import filters
from skimage import measure
from skimage import morphology
# from skimage import feature
from skimage import segmentation

# from scipy.ndimage import measurements
# from scipy import signal
from scipy import ndimage as ndi

visualize = True  # set true if you want to see the intermediate images
check = False  # set true if you want to confirm each crell manually
img_raw = io.imread('Photos/DRAQ50.1mMDAPI200ngul_Cruzella_TP13G-01.tif')  # the tif file to analyse

img_ch0 = np.max(img_raw[:, :, :, 1], axis=0)  # DAPI, Be fcng careful with channels, they change!
img_ch1 = np.max(img_raw[:, :, :, 0], axis=0)  # trans, rewrite only the number in []
img_ch2 = np.max(img_raw[:, :, :, 2], axis=0)  # Draq5

img_ch0_sum = np.sum(img_raw[:, :, :, 0], axis=0) / len(img_raw[:, :, :, 0])  # DAPI sum for measuring dna
img_ch2_sum = np.sum(img_raw[:, :, :, 2], axis=0) / len(img_raw[:, :, :, 2])  # Draq5

dna_mass = []

model = [{'sigma': 2, 'grad': 1, 'erosion': 10, 'dilation': 10, 'threshold': 1.5},
         {'sigma': 2, 'grad': 3, 'erosion': 10, 'dilation': 10, 'threshold': 1},
         {'sigma': 1, 'grad': 1, 'erosion': 9, 'dilation': 10, 'threshold': 1.5}]
model_num = 1  # use to adjust parameters TODO manual selection of a cell

if visualize:  # draw max values in 3 channels
    plt.figure(figsize=(30, 30))
    ax0 = plt.subplot(311)
    ax0.imshow(img_ch0, cmap='jet')
    ax0.set_title('DAPI')
    ax0.axis('off')

    ax1 = plt.subplot(312)
    ax1.imshow(img_ch1)
    ax1.set_title('trans')
    ax1.axis('off')

    ax2 = plt.subplot(313)
    ax2.imshow(img_ch2, cmap='jet')
    ax2.set_title('Draq5')
    ax2.axis('off')

    plt.show()

# Processing to contrast the cells, params set earlier
img_ch1_filt = filters.gaussian(img_ch1, sigma=model[model_num]['sigma'])
img_ch1_filt = img_ch1_filt / np.max(np.abs(img_ch1_filt))
img_ch1_filt = filters.rank.gradient(img_ch1_filt, morphology.disk(model[model_num]['grad']))

mask = img_ch1_filt > filters.threshold_otsu(img_ch1_filt)
mask = ndi.binary_fill_holes(mask)
mask = segmentation.clear_border(mask)

mask = morphology.erosion(mask, footprint=morphology.disk(model[model_num]['erosion']))
mask = morphology.dilation(mask, footprint=morphology.disk(model[model_num]['dilation']))

labels, labels_num = ndi.label(mask)

if visualize:
    plt.figure(figsize=(10, 10))
    plt.imshow(img_ch1, cmap='Greys')
    plt.imshow(ma.masked_where(~mask, labels), alpha=.3)

ctrl_fluo_img = img_ch0 + img_ch2
ctrl_fluo_mask = ctrl_fluo_img > filters.threshold_otsu(ctrl_fluo_img) * model[model_num]['threshold']

sums = ndi.sum_labels(ctrl_fluo_mask, labels, np.arange(labels_num + 1))

connected = sums > 0
debris_mask = connected[labels]

fin_mask = np.copy(mask)
fin_mask[~debris_mask] = 0

cells_labels, cells_num = ndi.label(fin_mask)

print(cells_num)

if visualize:
    plt.figure(figsize=(15, 15))
    plt.imshow(ctrl_fluo_mask, cmap='Greys')
    plt.imshow(ma.masked_where(~fin_mask, fin_mask), alpha=.3)
    plt.imshow(ma.masked_where(~debris_mask, debris_mask), cmap='jet', alpha=.3)


class OneCell:
    def __init__(self, dapi_img, draq_img, trans_img, cell_mask,
                 dapi_img_sum, draq_img_sum, visualize_bool, cell_area, cell_id):
        self.kn_label = None
        self.kn_mask = None
        self.n_mask = None
        self.k_mask = None
        self.k_mask_raw = None
        self.dapi_img = dapi_img
        self.draq_img = draq_img
        self.trans_img = trans_img
        self.cell_mask = cell_mask
        self.dapi_img_sum = dapi_img_sum
        self.draq_img_sum = draq_img_sum
        self.visualize = visualize_bool
        self.area = cell_area
        self.ID = cell_id
        self.dapi_norm = (self.dapi_img - np.min(self.dapi_img)) / (np.max(self.dapi_img) - np.min(self.dapi_img))
        self.draq_norm = (self.draq_img - np.min(self.draq_img)) / (np.max(self.draq_img) - np.min(self.draq_img))
        zeros = np.zeros_like(self.dapi_norm)

        self.rgb_overlay = np.stack([self.dapi_norm, self.draq_norm, zeros], axis=-1)  # red - DAPI, green - Draq5

    def show_imgs(self):

        plt.figure(figsize=(15, 15))

        ax0 = plt.subplot(221)
        ax0.imshow(self.trans_img, cmap='Greys')
        ax0.imshow(ma.masked_where(~self.cell_mask, self.cell_mask), cmap='jet', alpha=.3)
        ax0.set_title('Trans ch. with cell mask overlay')
        ax0.axis('off')

        ax1 = plt.subplot(222)
        ax1.imshow(self.rgb_overlay)
        ax1.set_title('ID=' + str(self.ID) + ', DAPI-red, Draq5-green')
        ax1.axis('off')

        plt.show()

    def create_kineto_mask(self):
        global dna_mass
        diff = self.dapi_norm - self.draq_norm

        self.k_mask_raw = diff > filters.threshold_otsu(diff)

        diff_masked = ma.masked_where(~self.k_mask_raw, diff)

        self.k_mask = diff > filters.threshold_otsu(diff_masked.compressed())

        self.n_mask = self.draq_norm > filters.threshold_otsu(self.draq_norm)

        # thresholds = filters.threshold_multiotsu(diff)
        # self.k_regions = np.digitize(diff, bins=thresholds)

        if np.sum(self.k_mask) > np.size(diff) * 0.5:
            print("Bad mask")

        nucleus = ma.masked_where(~self.n_mask, self.draq_img_sum)  # use the sum images, not max
        kinetoplast = ma.masked_where(~self.k_mask, self.draq_img_sum)
        nucleus_dna = ma.sum(nucleus)  # count the sum of intencities of the nucleus
        kineto_dna = ma.sum(kinetoplast)  # count the sum of intencities of the kinetoplast
        print(nucleus_dna, kineto_dna, kineto_dna / nucleus_dna)

        dna_mass.append((kineto_dna, nucleus_dna, kineto_dna / nucleus_dna, self.area, self.ID))
        # the array with the results

        if self.visualize:
            plt.figure(figsize=(10, 10))
            plt.imshow(nucleus, cmap='Blues')
            plt.imshow(kinetoplast, cmap='Greens')

            plt.title('ID=' + str(self.ID) + ', Blue-nucleus, Green - kineto')
            plt.axis('off')

            plt.show()


if visualize:
    fig, ax = plt.subplots(figsize=(15, 15))
    ax.imshow(img_ch1, cmap='Greys')

cells_list = []
area_mass = []
for region in measure.regionprops(cells_labels):
    minr, minc, maxr, maxc = region.bbox
    area = region.area
    area_mass.append((area, minr, minc, maxr, maxc))  # measure the area of the cells

area_mass = list(
    filter(
        lambda x: np.average(area_mass) - 5 * np.std(area_mass) < x[0] < np.average(area_mass) + 5 * np.std(area_mass),
        area_mass))  # to exclude anything too big or too small, TODO clasterise into putative cell types
if check:
    plt.ion()
ID = 0  # give an ID to each cell
for region in measure.regionprops(cells_labels):
    minr, minc, maxr, maxc = region.bbox
    if (region.area, minr, minc, maxr, maxc) not in area_mass:
        continue  # to exclude anything too big or too small
    if visualize:
        rect = mpatches.Rectangle((minc, minr), maxc - minc, maxr - minr, fill=False, edgecolor='red', linewidth=2)
        ax.add_patch(rect)  # doesn't matter, as visualize is the same condition, though ugly
        if check:
            plt.pause(2)  # let it think
            if input('Y/N:    ') == 'N':  # answer if it is a cell to count
                rect.remove()
                plt.ioff()
                continue
        # plt.ioff()
    ID += 1
    # print(ID)
    cells_list.append(OneCell(dapi_img=img_ch0[minr:maxr, minc:maxc],
                              draq_img=img_ch2[minr:maxr, minc:maxc],
                              trans_img=img_ch1[minr:maxr, minc:maxc],
                              cell_mask=fin_mask[minr:maxr, minc:maxc],
                              dapi_img_sum=img_ch0_sum[minr:maxr, minc:maxc],
                              draq_img_sum=img_ch2_sum[minr:maxr, minc:maxc],
                              visualize_bool=visualize,
                              cell_area=region.area, cell_id=ID))

print(len(cells_list))

if visualize:
    ax.set_axis_off()
    plt.tight_layout()
    if check:
        plt.ioff()
    plt.show()

for cell in cells_list:
    if visualize:
        cell.show_imgs()
    cell.create_kineto_mask()

dna_mass.sort(key=lambda i: i[0])  # /i[3] ** (3 / 2)
print(dna_mass)

# somehow plot the results
fig, ax = plt.subplots()
dots = [ax.annotate(str(i[4]) + ', ' + str(round(i[3] ** (3 / 2))), (i[0], i[1])) for i in dna_mass]
x = [i[0] for i in dna_mass]
y = [i[1] for i in dna_mass]
plt.plot(x, y, ':o')
plt.title('Nuclear DNA / Kinetoplast DNA')
# plt.plot(x, y)
plt.show()
