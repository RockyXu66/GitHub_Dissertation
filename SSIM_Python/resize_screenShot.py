# import numpy as np
# import cv2
# from matplotlib import pyplot as plt

#  # Load an color image in grayscale
# img = cv2.imread("test.jpg", 0) # 1-rgb, 0-grayscale, -1  Loads image as such including alpha channel
# img = cv2.resize(img, (0, 0), fx=0.5, fy=0.5)
# # cv2.namedWindow('image', cv2.WINDOW_NORMAL)
# # cv2.imshow('image',img)

# plt.imshow(img, cmap = 'gray', interpolation = 'bicubic')
# plt.xticks([]), plt.yticks([])  # to hide tick values on X and Y axis
# plt.show()


# # k = cv2.waitKey(0)
# # if k == 27: 		# wait for ESC key to exit
# # 	cv2.destroyAllWindows()
# # elif k == ord('s'):  # wait for 's' key to save and exit
# # 	cv2.imwrite('grayImage.jpg', img)
# # 	cv2.destroyAllWindows()

# import cv2
# for i in range(1, 37):
# 	print(i)
# 	img = cv2.imread("new_Model_screenShot/"+str(i)+".png")
# 	crop_img = img[140:160, 140:160] # Crop from x, y, w, h -> 100, 200, 300, 400
# # NOTE: its img[y: y + h, x: x + w] and *not* img[x: x + w, y: y + h]
# # cv2.imshow("cropped", crop_img)
# 	cv2.imwrite('cell_screenShot/'+str(i)+'.png', crop_img)
# 	# cv2.imshow(str(i), crop_img)	

# cv2.waitKey(0)

# import cv2

# testImage = cv2.imread("images/RGB46_127_179.png")
# # Opencv always uses BGR channel order
# p = testImage[0,299, 2]
# print p

# refImage = cv2.imread("new_Model_screenShot/1.png")
# cv2.imshow('reference image', refImage)

# file = open("testfile.txt","w") 
 
# file.write("Hello World\n") 
# file.write("This is our new text file") 
# file.write("and this is another line.") 
# file.write("Why? Because we can.") 
 
# file.close() 

# if cv2.waitKey(0) & 0xFF == ord('q'):
# 	cv2.destroyAllWindows()


import cv2

# for i in range(900):
# 	img = cv2.imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/head_900(resolution1080)/head"+str(i+1)+".png")
# 	size = (540,540)
# 	img = cv2.resize(img, size, interpolation=cv2.INTER_AREA)

# # cv2.imshow("resize", img)
# 	cv2.imwrite("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/head_900(resolution540)/head"+str(i+1)+".png", img)
# # cv2.waitKey(0)
# print("Finish resizing")


img = cv2.imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/ResultImages/head(720)_cell10_10/head11.png")
img = img[360:420, 360:420]
cv2.imwrite("cells10_10.png", img)
# cv2.imshow("cells", img)
# cv2.waitKey(0)




























