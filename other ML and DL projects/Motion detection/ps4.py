"""Problem Set 4: Motion Detection"""

import cv2
import numpy as np
import scipy as sp
from scipy import signal


# Utility function
def read_video(video_file, show=False):
    """Reads a video file and outputs a list of consecuative frames
  Args:
      image (string): Video file path
      show (bool):    Visualize the input video. WARNING doesn't work in
                      notebooks
  Returns:
      list(numpy.ndarray): list of frames
  """
    frames = []
    cap = cv2.VideoCapture(video_file)
    while (cap.isOpened()):
        ret, frame = cap.read()
        if not ret:
            break
        frames.append(frame)

        # Opens a new window and displays the input
        if show:
            cv2.imshow("input", frame)
            # Frames are read by intervals of 1 millisecond. The
            # programs breaks out of the while loop when the
            # user presses the 'q' key
            if cv2.waitKey(1) & 0xFF == ord('q'):
                break

    # The following frees up resources and
    # closes all windows
    cap.release()
    if show:
        cv2.destroyAllWindows()
    return frames
    
def normalize_and_scale(image_in, scale_range=(0, 255)):
    """Normalizes and scales an image to a given range [0, 255].

    Utility function. There is no need to modify it.

    Args:
        image_in (numpy.array): input image.
        scale_range (tuple): range values (min, max). Default set to [0, 255].

    Returns:
        numpy.array: output image.
    """
    image_out = np.zeros(image_in.shape)
    cv2.normalize(image_in, image_out, alpha=scale_range[0],
                  beta=scale_range[1], norm_type=cv2.NORM_MINMAX)

    return image_out


# Assignment code
def gradient_x(image):
    """Computes image gradient in X direction.

    Use cv2.Sobel to help you with this function. Additionally you
    should set cv2.Sobel's 'scale' parameter to one eighth and ksize
    to 3.

    Args:
        image (numpy.array): grayscale floating-point image with values in [0.0, 1.0].

    Returns:
        numpy.array: image gradient in the X direction. Output
                     from cv2.Sobel.
    """
    imgradx = cv2.Sobel(image, cv2.CV_64F,1,0,ksize=3, scale = 1/8)
    return imgradx
    raise NotImplementedError


def gradient_y(image):
    """Computes image gradient in Y direction.
    Use cv2.Sobel to help you with this function. Additionally you
    should set cv2.Sobel's 'scale' parameter to one eighth and ksize
    to 3.

    Args:
        image (numpy.array): grayscale floating-point image with values in [0.0, 1.0].

    Returns:
        numpy.array: image gradient in the Y direction.
                     Output from cv2.Sobel.
    """
    imgrady = cv2.Sobel(image, cv2.CV_64F,0,1,ksize=3, scale = 1/8)
    return imgrady
    raise NotImplementedError


def optic_flow_lk(img_a, img_b, k_size, k_type, sigma=1):
    """Computes optic flow using the Lucas-Kanade method.

    For efficiency, you should apply a convolution-based method.
    Note: Implement this method using the instructions in the lectures
    and the documentation.
    You are not allowed to use any OpenCV functions that are related
    to Optic Flow.
    Args:
        img_a (numpy.array): grayscale floating-point image with
                             values in [0.0, 1.0].
        img_b (numpy.array): grayscale floating-point image with
                             values in [0.0, 1.0].
        k_size (int): size of averaging kernel to use for weighted
                      averages. Here we assume the kernel window is a
                      square so you will use the same value for both
                      width and height.
        k_type (str): type of kernel to use for weighted averaging,
                      'uniform' or 'gaussian'. By uniform we mean a
                      kernel with the only ones divided by k_size**2.
                      To implement a Gaussian kernel use
                      cv2.getGaussianKernel. The autograder will use
                      'uniform'.
        sigma (float): sigma value if gaussian is chosen. Default
                       value set to 1 because the autograder does not
                       use this parameter.
    Returns:
        tuple: 2-element tuple containing:
            U (numpy.array): raw displacement (in pixels) along
                             X-axis, same size as the input images,
                             floating-point type.
            V (numpy.array): raw displacement (in pixels) along
                             Y-axis, same size and type as U.
    """
    Ix = gradient_x(img_a)
    Iy = gradient_y(img_a)
    if k_type == 'gaussian':
        k = cv2.getGaussianKernel(ksize, sigma)
        gk = k * k.transpose(1, 0)
    if k_type == 'uniform':
        gk = np.ones((k_size,k_size))
    It = img_b - img_a
    xx = Ix*Ix
    xy = Ix*Iy
    yy = Iy*Iy
    xt = -Ix*It
    yt = -Iy*It
    a = sp.signal.convolve2d(xx, gk, mode='same', boundary='fill')
    b = sp.signal.convolve2d(xy, gk, mode='same', boundary='fill')
    c = b
    d = sp.signal.convolve2d(yy, gk, mode='same', boundary='fill')
    e = sp.signal.convolve2d(xt, gk, mode='same', boundary='fill')
    f = sp.signal.convolve2d(yt, gk, mode='same', boundary='fill')
    print(a.shape)
    u = img_a*0
    v = img_a*0
    for i in range(0,img_a.shape[0]):
        for j in range(0,img_a.shape[1]):
            res = np.matmul(np.linalg.pinv(np.array([[a[i][j],b[i][j]],[c[i][j],d[i][j]]])),np.array([[e[i][j]],[f[i][j]]]))
            u[i][j] = res[0]
            v[i][j] = res[1]
    return u,v
    raise NotImplementedError


def reduce_image(image):
    """Reduces an image to half its shape.

    The autograder will pass images with even width and height. It is
    up to you to determine values with odd dimensions. For example the
    output image can be the result of rounding up the division by 2:
    (13, 19) -> (7, 10)

    For simplicity and efficiency, implement a convolution-based
    method using the 5-tap separable filter.

    Follow the process shown in the lecture 6B-L3. Also refer to:
    -  Burt, P. J., and Adelson, E. H. (1983). The Laplacian Pyramid
       as a Compact Image Code
    You can find the link in the problem set instructions.

    Args:
        image (numpy.array): grayscale floating-point image, values in
                             [0.0, 1.0].

    Returns:
        numpy.array: output image with half the shape, same type as the
                     input image.
    """
    inter = image.copy()
    out = image.copy()
    hf = np.array([1/16, 1/4, 3/8, 1/4, 1/16])
    vf = np.array([[1/16], [1/4], [3/8], [1/4], [1/16]])
    # for i in range(0,image.shape[0]):
    #     inter[i,:]=sp.signal.convolve(image[i,:],hf,mode='same')
    # for j in range(0,image.shape[1]):
    #     col = image[:,j]
    #     col = col.reshape((image.shape[0],1))
    #     inter[:,[j]]=sp.signal.convolve(col,vf,mode='same')
    inter = cv2.sepFilter2D(inter, -1, hf, vf, borderType=cv2.BORDER_REFLECT_101)
    for m in range(0,int(np.round((image.shape[0]+0.01)/2))):
        for n in range(0,int(np.round((image.shape[1]+0.01)/2))):
            out[m,n] = inter[2*m,2*n]
    return out[0:int(np.round((image.shape[0]+0.01)/2)),0:int(np.round((image.shape[1]+0.01)/2))]
    raise NotImplementedError


def gaussian_pyramid(image, levels):
    """Creates a Gaussian pyramid of a given image.

    This method uses reduce_image() at each level. Each image is
    stored in a list of length equal the number of levels.

    The first element in the list ([0]) should contain the input
    image. All other levels contain a reduced version of the previous
    level.

    All images in the pyramid should floating-point with values in

    Args:
        image (numpy.array): grayscale floating-point image, values
                             in [0.0, 1.0].
        levels (int): number of levels in the resulting pyramid.

    Returns:
        list: Gaussian pyramid, list of numpy.arrays.
    """
    imglist = []
    im = image.copy()
    for l in range(0,levels):
        imglist.append(im)
        im = reduce_image(im)
    return imglist
    raise NotImplementedError


def create_combined_img(img_list):
    """Stacks images from the input pyramid list side-by-side.

    Ordering should be large to small from left to right.

    See the problem set instructions for a reference on how the output
    should look like.

    Make sure you call normalize_and_scale() for each image in the
    pyramid when populating img_out.

    Args:
        img_list (list): list with pyramid images.

    Returns:
        numpy.array: output image with the pyramid images stacked
                     from left to right.
    """
    l = np.shape(img_list)
    l = l[0]
    h = img_list[0].shape[0]
    w = img_list[0].shape[1]
    out = np.zeros((h,w*2))
    rs = 0
    cs = 0
    re = 0
    ce = 0
    for i in range(0,l):
        pic = normalize_and_scale(img_list[i])
        ce = ce + pic.shape[1]
        re = pic.shape[0]
        out[rs:re,cs:ce] = pic
        cs = cs + pic.shape[1]
        print(cs,ce)
    return out[:,0:ce]
    raise NotImplementedError


def expand_image(image):
    """Expands an image doubling its width and height.

    For simplicity and efficiency, implement a convolution-based
    method using the 5-tap separable filter.

    Follow the process shown in the lecture 6B-L3. Also refer to:
    -  Burt, P. J., and Adelson, E. H. (1983). The Laplacian Pyramid
       as a Compact Image Code

    You can find the link in the problem set instructions.

    Args:
        image (numpy.array): grayscale floating-point image, values
                             in [0.0, 1.0].

    Returns:
        numpy.array: same type as 'image' with the doubled height and
                     width.
    """
    h = image.shape[0]
    w = image.shape[1]
    pout = np.zeros((h*2,w*2))
    inter = image.copy()
    for m in range(0,int(np.round((pout.shape[0]+0.01)/2))):
        for n in range(0,int(np.round((pout.shape[1]+0.01)/2))):
            pout[m*2,n*2] = inter[m,n]
    hf = np.array([1/8, 1/2, 6/8, 1/2, 1/8])
    vf = np.array([[1/8], [1/2], [6/8], [1/2], [1/8]])
    out = cv2.sepFilter2D(pout, -1, hf, vf, borderType=cv2.BORDER_REFLECT_101)
    return out
    raise NotImplementedError


def laplacian_pyramid(g_pyr):
    """Creates a Laplacian pyramid from a given Gaussian pyramid.

    This method uses expand_image() at each level.

    Args:
        g_pyr (list): Gaussian pyramid, returned by gaussian_pyramid().

    Returns:
        list: Laplacian pyramid, with l_pyr[-1] = g_pyr[-1].
    """
    l = np.shape(g_pyr)
    levels = l[0]
    imglist = []
    im = g_pyr[-1]
    imglist.insert(0,im)
    for l in range(0,levels-1):
        gsubtr1 = g_pyr[-1-l]
        im = expand_image(gsubtr1)
        gsubtr2 = g_pyr[-2-l]
        im_out = gsubtr2 - im[0:gsubtr2.shape[0],0:gsubtr2.shape[1]]
#         im_out = gsubtr - im
        imglist.insert(0,im_out)
    return imglist
    raise NotImplementedError


def warp(image, U, V, interpolation, border_mode):
    """Warps image using X and Y displacements (U and V).

    This function uses cv2.remap. The autograder will use cubic
    interpolation and the BORDER_REFLECT101 border mode. You may
    change this to work with the problem set images.

    See the cv2.remap documentation to read more about border and
    interpolation methods.

    Args:
        image (numpy.array): grayscale floating-point image, values
                             in [0.0, 1.0].
        U (numpy.array): displacement (in pixels) along X-axis.
        V (numpy.array): displacement (in pixels) along Y-axis.
        interpolation (Inter): interpolation method used in cv2.remap.
        border_mode (BorderType): pixel extrapolation method used in
                                  cv2.remap.

    Returns:
        numpy.array: warped image, such that
                     warped[y, x] = image[y + V[y, x], x + U[y, x]]
    """
    M, N = image.shape
    X, Y = np.meshgrid(range(N), range(M))
    X_map = np.zeros((image.shape[0], image.shape[1]), dtype=np.float32)
    Y_map = np.zeros((image.shape[0], image.shape[1]), dtype=np.float32)
    X_map = X + U
    Y_map = Y + V
    src = image.astype(np.float32)
    X_map = X_map.astype(np.float32)
    Y_map = Y_map.astype(np.float32)
    dst = cv2.remap(src, X_map, Y_map,interpolation=interpolation,borderMode=border_mode)
    return dst
    raise NotImplementedError


def hierarchical_lk(img_a, img_b, levels, k_size, k_type, sigma, interpolation,
                    border_mode):
    """Computes the optic flow using Hierarchical Lucas-Kanade.

    This method should use reduce_image(), expand_image(), warp(),
    and optic_flow_lk().

    Args:
        img_a (numpy.array): grayscale floating-point image, values in
                             [0.0, 1.0].
        img_b (numpy.array): grayscale floating-point image, values in
                             [0.0, 1.0].
        levels (int): Number of levels.
        k_size (int): parameter to be passed to optic_flow_lk.
        k_type (str): parameter to be passed to optic_flow_lk.
        sigma (float): parameter to be passed to optic_flow_lk.
        interpolation (Inter): parameter to be passed to warp.
        border_mode (BorderType): parameter to be passed to warp.

    Returns:
        tuple: 2-element tuple containing:
            U (numpy.array): raw displacement (in pixels) along X-axis,
                             same size as the input images,
                             floating-point type.
            V (numpy.array): raw displacement (in pixels) along Y-axis,
                             same size and type as U.
    """
    pa=gaussian_pyramid(img_a, levels)
    pb=gaussian_pyramid(img_b, levels)
    la=pa[-1]
    lb=pb[-1]
    print(pa[0].shape)
    ut,vt = optic_flow_lk(lb, la, k_size, k_type)
    ut = 2*expand_image(ut)
    vt = 2*expand_image(vt)
    for i in range(levels-2,0,-1):
        warped = warp(pb[i], -ut, -vt, interpolation, border_mode)
        print('warped',warped.shape)
        u,v = optic_flow_lk(warped, pa[i], k_size, k_type)

        print('u',u.shape)
        print('ut',ut.shape)
        ut = ut + u
        vt = vt + v

        ut = 2*expand_image(ut)
        vt = 2*expand_image(vt)
#         if (ut.shape != pb[i-1].shape or vt.shape != pb[i-1].shape) and i!=0:
#             ut = ut[0:pb[i-1].shape[0], 0:pb[i-1].shape[1]]
#             vt = vt[0:pb[i-1].shape[0], 0:pb[i-1].shape[1]]

    warped = warp(pb[0], -ut, -vt, interpolation, border_mode)
    u,v = optic_flow_lk(warped, pa[0], k_size, k_type)
    ut = ut + u
    vt = vt + v
    return -ut,-vt
    raise NotImplementedError

def classify_video(images):


    return
