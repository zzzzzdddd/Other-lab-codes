"""
2d Max Pooling Module.  (c) 2021 Georgia Tech

Copyright 2021, Georgia Institute of Technology (Georgia Tech)
Atlanta, Georgia 30332
All Rights Reserved

Template code for CS 7643 Deep Learning

Georgia Tech asserts copyright ownership of this template and all derivative
works, including solutions to the projects assigned in this course. Students
and other users of this template code are advised not to share it with others
or to make it available on publicly viewable websites including repositories
such as Github, Bitbucket, and Gitlab.  This copyright statement should
not be removed or edited.

Sharing solutions with current or future students of CS 7643 Deep Learning is
prohibited and subject to being investigated as a GT honor code violation.

-----do not edit anything above this line---
"""

import numpy as np


class MaxPooling:
    """
    Max Pooling of input
    """

    def __init__(self, kernel_size, stride):
        self.kernel_size = kernel_size
        self.stride = stride
        self.cache = None
        self.dx = None

    def forward(self, x):
        """
        Forward pass of max pooling
        :param x: input, (N, C, H, W)
        :return: The output by max pooling with kernel_size and stride
        """
        out = None
        #############################################################################
        # TODO: Implement the max pooling forward pass.                             #
        # Hint:                                                                     #
        #       1) You may implement the process with loops                         #
        #############################################################################
        pic = x.copy()
        # print(x.shape, pic.shape)
        hp = int((pic.shape[2] - self.kernel_size)/self.stride) + 1
        wp = int((pic.shape[3] - self.kernel_size)/self.stride) + 1
        out = np.zeros([pic.shape[0], pic.shape[1], hp, wp])
        # for oc in range(pic.shape[1]):
        for i in range(hp):
            for j in range(wp):
                out[:,:,i,j] = np.amax(pic[:,:,i*self.stride:i*self.stride+self.kernel_size,j*self.stride:j*self.stride+self.kernel_size], axis=(2,3))
        H_out = hp
        W_out = wp
        #############################################################################
        #                              END OF YOUR CODE                             #
        #############################################################################
        self.cache = (x, H_out, W_out)
        return out

    def backward(self, dout):
        """
        Backward pass of max pooling
        :param dout: Upstream derivatives
        :return:
        """
        x, H_out, W_out = self.cache
        #############################################################################
        # TODO: Implement the max pooling backward pass.                            #

        #############################################################################
        self.dx = np.zeros(x.shape)
        for n in range(x.shape[0]):
            for c in range(x.shape[1]):
                for i in range(H_out):
                    for j in range(W_out):
                        maxi, maxj = np.where(np.amax(x[n,c,i*self.stride:i*self.stride+self.kernel_size,j*self.stride:j*self.stride+self.kernel_size])\
                                              ==x[n,c,i*self.stride:i*self.stride+self.kernel_size,j*self.stride:j*self.stride+self.kernel_size])
                        self.dx[n,c,i*self.stride:i*self.stride+self.kernel_size,j*self.stride:j*self.stride+self.kernel_size][maxi[0], maxj[0]] = dout[n,c,i,j]
        #############################################################################
        #                              END OF YOUR CODE                             #
        #############################################################################
