"""
LSTM model.  (c) 2021 Georgia Tech

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
import torch
import torch.nn as nn


class LSTM(nn.Module):
    # An implementation of naive LSTM using Pytorch Linear layers and activations
    # You will need to complete the class init function, forward function and weight initialization

    def __init__(self, input_size, hidden_size):
        """ Init function for LSTM class
            Args:
                input_size (int): the number of features in the inputs.
                hidden_size (int): the size of the hidden layer
            Returns: 
                None
        """
        super(LSTM, self).__init__()

        self.input_size = input_size
        self.hidden_size = hidden_size

        ################################################################################
        # TODO:                                                                        #
        #   Declare LSTM weights and attributes in order specified below to pass GS.   #
        #   You should include weights and biases regarding using nn.Parameter:        #
        #       1) i_t: input gate                                                     #
        #       2) f_t: forget gate                                                    #
        #       3) g_t: cell gate, or the tilded cell state                            #
        #       4) o_t: output gate                                                    #
        #   for each equation above, initialize the weights,biases for input prior     #
        #   to weights, biases for hidden.                                             #
        #   You also need to include correct activation functions                      #
        ################################################################################

        # i_t: input gate
        self.i_xweights = nn.Parameter(torch.randn(self.input_size, self.hidden_size))
        self.i_xbias = nn.Parameter(torch.zeros(self.hidden_size))
        self.i_hweights = nn.Parameter(torch.randn(self.hidden_size, self.hidden_size))
        self.i_hbias = nn.Parameter(torch.zeros(self.hidden_size))
        # f_t: the forget gate
        self.f_xweights = nn.Parameter(torch.randn(self.input_size, self.hidden_size))
        self.f_xbias = nn.Parameter(torch.zeros(self.hidden_size))
        self.f_hweights = nn.Parameter(torch.randn(self.hidden_size, self.hidden_size))
        self.f_hbias = nn.Parameter(torch.zeros(self.hidden_size))
        # g_t: the cell gate
        self.g_xweights = nn.Parameter(torch.randn(self.input_size, self.hidden_size))
        self.g_xbias = nn.Parameter(torch.zeros(self.hidden_size))
        self.g_hweights = nn.Parameter(torch.randn(self.hidden_size, self.hidden_size))
        self.g_hbias = nn.Parameter(torch.zeros(self.hidden_size))
        # o_t: the output gate
        self.o_xweights = nn.Parameter(torch.randn(self.input_size, self.hidden_size))
        self.o_xbias = nn.Parameter(torch.zeros(self.hidden_size))
        self.o_hweights = nn.Parameter(torch.randn(self.hidden_size, self.hidden_size))
        self.o_hbias = nn.Parameter(torch.zeros(self.hidden_size))

        ################################################################################
        #                              END OF YOUR CODE                                #
        ################################################################################
        self.init_hidden()

    def init_hidden(self):
        for p in self.parameters():
            if p.data.ndimension() >= 2:
                nn.init.xavier_uniform_(p.data)
            else:
                nn.init.zeros_(p.data)

    def forward(self, x: torch.Tensor, init_states=None):
        """Assumes x is of shape (batch, sequence, feature)"""

        ################################################################################
        # TODO:                                                                        #
        #   Implement the forward pass of LSTM. Please refer to the equations in the   #
        #   corresponding section of jupyter notebook. Iterate through all the time    #
        #   steps and return only the hidden and cell state, h_t and c_t.              # 
        #   Note that this time you are also iterating over all of the time steps.     #
        ################################################################################
        h_t, c_t = None, None
        # print(self.input_size, self.hidden_size)
        # print(x.size())
        ht_1 = torch.zeros(x.size(0), self.hidden_size)
        ct_1 = torch.zeros(x.size(0), self.hidden_size)

        sig = nn.Sigmoid()
        tanh = nn.Tanh()
        for i in range(x.size(1)):
            xt = x[:,i,:].T
            # print(xt)
            # print(self.i_xweights.data.size())
            # print(xt.size())
            # print(self.i_xbias.data.size())
            # print(sig(torch.sum(self.i_xweights.data.T@xt + self.i_xbias.data + self.i_hweights.data.T@ht_1 + self.i_hbias.data)))
            it = sig(self.i_xweights.data.T@xt + self.i_xbias.data + self.i_hweights.data.T@ht_1 + self.i_hbias.data)
            ft = sig(self.f_xweights.data.T@xt + self.f_xbias.data + self.f_hweights.data.T@ht_1 + self.f_hbias.data)
            gt = tanh(self.g_xweights.data.T@xt + self.g_xbias.data + self.g_hweights.data.T@ht_1 + self.g_hbias.data)
            ot = sig(self.o_xweights.data.T@xt + self.o_xbias.data + self.o_hweights.data.T@ht_1 + self.o_hbias.data)

            c_t = ft * ct_1 + it * gt
            h_t = ot * tanh(c_t)
            ht_1 = h_t
            ct_1 = c_t
        ################################################################################
        #                              END OF YOUR CODE                                #
        ################################################################################
        # print(h_t, c_t)
        return (h_t.T, c_t.T)
