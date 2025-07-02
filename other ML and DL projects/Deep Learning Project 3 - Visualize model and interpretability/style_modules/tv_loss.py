import torch
import torch.nn as nn

class TotalVariationLoss(nn.Module):
    def forward(self, img, tv_weight):
        """
            Compute total variation loss.

            Inputs:
            - img: PyTorch Variable of shape (1, 3, H, W) holding an input image.
            - tv_weight: Scalar giving the weight w_t to use for the TV loss.

            Returns:
            - loss: PyTorch Variable holding a scalar giving the total variation loss
              for img weighted by tv_weight.
            """

        ##############################################################################
        # TODO: Implement total varation loss function                               #
        # Please pay attention to use torch tensor math function to finish it.       #
        # Otherwise, you may run into the issues later that dynamic graph is broken  #
        # and gradient can not be derived.                                           #
        ##############################################################################
        H = list(img.shape)[2]
        W = list(img.shape)[3]
        sig1 = torch.sum(torch.square(img[:,:, :, 0:W-1] - img[:,:,:, 1:W]), (2,3))
        sig2 = torch.sum(torch.square(img[:,:, 0:H-1, :] - img[:,:, 1:H, :]), (2,3))
        sig = torch.sum(sig1+sig2,1)
        loss = tv_weight * sig
        return loss
        ##############################################################################
        #                             END OF YOUR CODE                               #
        ##############################################################################
