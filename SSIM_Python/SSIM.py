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

# load the images -- the original, the original + contrast,
# and the original + photoshop
#original = cv2.imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/artifix_120/artifix11.png")
#contrast = cv2.imread("/Users/yinghanxu/Study/GitHub_Dissertation/PCA_opencv_cells/opencv tutorial/ResultImages/artifix11_cell40_40.png")

#original = cv2.imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/head_900(resolution1080)/head11.png");
#contrast = cv2.imread("/Users/yinghanxu/Study/GitHub_Dissertation/PCA_opencv_cells/opencv tutorial/ResultImages/head11_cell20_30.png")
original = cv2.imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/head_900(resolution720)/head11.png");
contrast = cv2.imread("/Users/yinghanxu/Study/GitHub_Dissertation/PCA_opencv_cells/opencv tutorial/ResultImages/head11(720)_cell20_25.png")

# convert the images to grayscale
original = cv2.cvtColor(original, cv2.COLOR_BGR2GRAY)
contrast = cv2.cvtColor(contrast, cv2.COLOR_BGR2GRAY)

## initialize the figure
#fig = plt.figure("Images")
images = ("Original", original), ("Contrast", contrast)

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
s = ssim(original, contrast, multichannel=True)
# print(s)

num_components = 20

cell_dimension = 20
image_width = 1080

original = cv2.imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/head_900(resolution"+str(image_width)+")/head11.png")
contrast = cv2.imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/ResultImages/head("+str(image_width)+")_cell"+str(cell_dimension)+"_"+str(num_components)+"/head11.png")
s = ssim(original, contrast, multichannel=True)
print s

# read_file = open("/Users/yinghanxu/Study/Dissertation_ResultData/SSIM/head"+str(image_width)+"_cells"+str(cell_dimension)+"_"+str(num_components)+".txt", "r")
# content = read_file.readline()
# t = float(content)
# read_file.close()

total = 0
for imageIndex in range(2,901):
    original = cv2.imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/head_900(resolution"+str(image_width)+")/head"+str(imageIndex)+".png")
    contrast = cv2.imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/ResultImages/head("+str(image_width)+")_cell"+str(cell_dimension)+"_"+str(num_components)+"/head"+str(imageIndex)+".png")
    s = ssim(original, contrast, multichannel=True)
    total += s
    print imageIndex
average = total/899.0
print average

file = open("/Users/yinghanxu/Study/GitHub_Dissertation/Experiments/SSIM/head"+str(image_width)+"_cells"+str(cell_dimension)+"_"+str(num_components)+".txt", "w")
file.write(str(average))
file.close()

print("Save average SSIM for 899 images. head"+str(image_width)+"_cells"+str(cell_dimension)+"_"+str(num_components))



















