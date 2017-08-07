from skimage.measure import compare_ssim as ssim
import matplotlib.pyplot as plt
import numpy as np
import cv2

def mse(imageA, imageB):
	# the 'Mean Squared Error' between the two images is the
	# sum of the squared difference between the two images;
	# NOTE: the two images must have the same dimension
	err = np.sum((imageA.astype("float") - imageB.astype("float")) ** 2)
	err /= float(imageA.shape[0] * imageA.shape[1])
	
	# return the MSE, the lower the error, the more "similar"
	# the two images are
	return err

def compare_images(imageA, imageB, title):
	# compute the mean squared error and structural similarity
	# index for the images
#	m = mse(imageA, imageB)
	s = ssim(imageA, imageB)

	# setup the figure
	fig = plt.figure(title)
	plt.suptitle("SSIM: %.4f" % (s))

	# show first image
	ax = fig.add_subplot(1, 2, 1)
	plt.imshow(imageA, cmap = plt.cm.gray)
	plt.axis("off")

	# show the second image
	ax = fig.add_subplot(1, 2, 2)
	plt.imshow(imageB, cmap = plt.cm.gray)
	plt.axis("off")

	# show the images
	plt.show()

def cal_SSIM(num_components, isSmoothed, cell_dimension, image_width):
	total = 0
	for imageIndex in range(2,901):
   		original = cv2.imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/head_900(resolution"+str(image_width)+")/head"+str(imageIndex)+".png")
   		if isSmoothed:
			contrast = cv2.imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/ResultImages_Smoothed/head("+str(image_width)+")_cell"+str(cell_dimension)+"_"+str(num_components)+"/head"+str(imageIndex)+".png")
   		else: 
   			contrast = cv2.imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/ResultImages/head("+str(image_width)+")_cell"+str(cell_dimension)+"_"+str(num_components)+"/head"+str(imageIndex)+".png")
   		s = ssim(original, contrast, multichannel=True)
   		total += s
   		print imageIndex, s

	average = total/899.0
	print average

	if isSmoothed:
		file = open("/Users/yinghanxu/Study/GitHub_Dissertation/Experiments/SSIM/head"+str(image_width)+"_cells"+str(cell_dimension)+"_"+str(num_components)+"_Smoothed.txt", "w")
	else:	
		file = open("/Users/yinghanxu/Study/GitHub_Dissertation/Experiments/SSIM/head"+str(image_width)+"_cells"+str(cell_dimension)+"_"+str(num_components)+".txt", "w")
	file.write(str(average))
	file.close()

	print("Save average SSIM for 899 images. head"+str(image_width)+"_cells"+str(cell_dimension)+"_"+str(num_components))


# load the images -- the original, the original + contrast,
# and the original + photoshop
#original = cv2.imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/artifix_120/artifix11.png")
#contrast = cv2.imread("/Users/yinghanxu/Study/GitHub_Dissertation/PCA_opencv_cells/opencv tutorial/ResultImages/artifix11_cell40_40.png")

#original = cv2.imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/head_900(resolution1080)/head11.png");
#contrast = cv2.imread("/Users/yinghanxu/Study/GitHub_Dissertation/PCA_opencv_cells/opencv tutorial/ResultImages/head11_cell20_30.png")
#original = cv2.imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/head_900(resolution720)/head11.png");
#contrast = cv2.imread("/Users/yinghanxu/Study/GitHub_Dissertation/PCA_opencv_cells/opencv tutorial/ResultImages/head11(720)_cell10_10.png")

# convert the images to grayscale
#original = cv2.cvtColor(original, cv2.COLOR_BGR2GRAY)
#contrast = cv2.cvtColor(contrast, cv2.COLOR_BGR2GRAY)

## initialize the figure
#fig = plt.figure("Images")
#images = ("Original", original), ("Contrast", contrast)

## loop over the images
#for (i, (name, image)) in enumerate(images):
#	# show the image
#	ax = fig.add_subplot(1, 3, i + 1)
#	ax.set_title(name)
#	plt.imshow(image, cmap = plt.cm.gray)
#	plt.axis("off")
#
## show the figure
#plt.show()

# compare the images
#compare_images(original, original, "Original vs. Original")
#compare_images(original, contrast, "Original vs. Contrast")
#s = ssim(original, contrast, multichannel=True)
#print(s)

# original = cv2.imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/head_900(resolution540)/head11.png");
# contrast = cv2.imread("/Users/yinghanxu/Study/GitHub_Dissertation/PCA_OpenCV_Prepocess/opencv tutorial/ResultImages/head11(540)_cell20_60.png");
# test = ssim(original, contrast, multichannel=True)
# print test

num_components = 4
isSmoothed = True
cell_dimension = 5
image_width = 720

# original = cv2.imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/head_900(resolution"+str(image_width)+")/head11.png")
# contrast = cv2.imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/ResultImages_Smoothed/head11_smoothedImage_blur_3.png");
# # contrast = cv2.imread("/Users/yinghanxu/Study/GitHub_Dissertation/PCA_OpenCV_Prepocess/opencv tutorial/ResultImages/head11("+str(image_width)+")_cell"+str(cell_dimension)+"_"+str(num_components)+".png");
# s = ssim(original, contrast, multichannel=True)
# print s

cal_SSIM(8, isSmoothed, cell_dimension, image_width)
cal_SSIM(7, isSmoothed, cell_dimension, image_width)
cal_SSIM(6, isSmoothed, cell_dimension, image_width)
cal_SSIM(5, isSmoothed, cell_dimension, image_width)
cal_SSIM(4, isSmoothed, cell_dimension, image_width)

# total = 0
# for imageIndex in range(2,901):
#    original = cv2.imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/head_900(resolution"+str(image_width)+")/head"+str(imageIndex)+".png")
#    if isSmoothed:
# 	contrast = cv2.imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/ResultImages_Smoothed/head("+str(image_width)+")_cell"+str(cell_dimension)+"_"+str(num_components)+"/head"+str(imageIndex)+".png")
#    else: 
#    	contrast = cv2.imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/ResultImages/head("+str(image_width)+")_cell"+str(cell_dimension)+"_"+str(num_components)+"/head"+str(imageIndex)+".png")
#    s = ssim(original, contrast, multichannel=True)
#    total += s
#    print imageIndex, s
# average = total/899.0
# print average

# if isSmoothed:
# 	file = open("/Users/yinghanxu/Study/GitHub_Dissertation/Experiments/SSIM/head"+str(image_width)+"_cells"+str(cell_dimension)+"_"+str(num_components)+"_Smoothed.txt", "w")
# else:	
# 	file = open("/Users/yinghanxu/Study/GitHub_Dissertation/Experiments/SSIM/head"+str(image_width)+"_cells"+str(cell_dimension)+"_"+str(num_components)+".txt", "w")
# file.write(str(average))
# file.close()

# print("Save average SSIM for 899 images. head"+str(image_width)+"_cells"+str(cell_dimension)+"_"+str(num_components))

###################
# read_file = open("/Users/yinghanxu/Study/Dissertation_ResultData/SSIM/head"+str(image_width)+"_cells"+str(cell_dimension)+"_"+str(num_components)+".txt", "r")
# content = read_file.readline()
# t = float(content)
# read_file.close()
















