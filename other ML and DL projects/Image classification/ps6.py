"""Problem Set 6: PCA, Boosting, Haar Features, Viola-Jones."""
import numpy as np
import cv2
import os

from helper_classes import WeakClassifier, VJ_Classifier


# assignment code
def load_images(folder, size=(32, 32)):
    """Load images to workspace.

    Args:
        folder (String): path to folder with images.
        size   (tuple): new image sizes

    Returns:
        tuple: two-element tuple containing:
            X (numpy.array): data matrix of flatten images
                             (row:observations, col:features) (float).
            y (numpy.array): 1D array of labels (int).
    """
    images_files = [f for f in os.listdir(folder)]
    imgs = [np.array(cv2.imread(os.path.join(folder, f), 0)) for f in images_files]
    imgs = [cv2.resize(x, size) for x in imgs]

    y = np.zeros(np.shape(images_files)[0])
    y.astype(int)
    for i in range(0,np.shape(images_files)[0]):
        imgs[i] = imgs[i].flatten()
        y[i] = int(images_files[i][7:9])
    X = np.array(imgs)
    return X, y



def split_dataset(X, y, p):
    """Split dataset into training and test sets.

    Let M be the number of images in X, select N random images that will
    compose the training data (see np.random.permutation). The images that
    were not selected (M - N) will be part of the test data. Record the labels
    accordingly.

    Args:
        X (numpy.array): 2D dataset.
        y (numpy.array): 1D array of labels (int).
        p (float): Decimal value that determines the percentage of the data
                   that will be the training data.

    Returns:
        tuple: Four-element tuple containing:
            Xtrain (numpy.array): Training data 2D array.
            ytrain (numpy.array): Training data labels.
            Xtest (numpy.array): Test data test 2D array.
            ytest (numpy.array): Test data labels.
    """
    np.random.seed(1)
    tst = np.random.choice(range(y.shape[0]), int(y.shape[0]*p), replace=False)
    i = np.array([False]*y.shape[0])
    i[tst] = True
    Xtrain = X[i]
    ytrain = y[i]
    Xtest = X[~i]
    ytest = y[~i]
    return (Xtrain, ytrain, Xtest, ytest)


def get_mean_face(x):
    """Return the mean face.

    Calculate the mean for each column.

    Args:
        x (numpy.array): array of flattened images.

    Returns:
        numpy.array: Mean face.
    """
    mf = np.zeros(x.shape[1])
    for i in range(0, mf.shape[0]):
        mf[i] = np.mean(x[:,i])
    return mf


def pca(X, k):
    """PCA Reduction method.

    Return the top k eigenvectors and eigenvalues using the covariance array
    obtained from X.


    Args:
        X (numpy.array): 2D data array of flatten images (row:observations,
                         col:features) (float).
        k (int): new dimension space

    Returns:
        tuple: two-element tuple containing
            eigenvectors (numpy.array): 2D array with the top k eigenvectors.
            eigenvalues (numpy.array): array with the top k eigenvalues.
    """
    # mf = get_mean_face(X)
    # A = X
    # sig = np.zeros((X.shape[1],X.shape[1]))
    # for i in range(0,X.shape[0]):
    #     A[i,:] = A[i,:] - mf

    # ATA = np.dot(A, np.transpose(A))
    # eigstuff = np.linalg.eigh(ATA)
    # vals = eigstuff[0]
    # vecs = eigstuff[1]
    # C = np.matmul(np.transpose(A),vecs)
    # D = np.matmul(vecs,ATA)
    # vals = vals[::-1]
    # vecs = C[::-1]
    # result = (np.transpose(vecs[0:k]),vals[0:k])

    A = X - np.array(X.mean(0), ndmin=2)
    C = np.matmul(np.transpose(A),A)
    vals, vecs = np.linalg.eigh(C)
    vals = vals[::-1][0:k]
    vecs = vecs.T[::-1][0:k]
    return (np.transpose(vecs), vals)

class Boosting:
    """Boosting classifier.

    Args:
        X (numpy.array): Data array of flattened images
                         (row:observations, col:features) (float).
        y (numpy.array): Labels array of shape (observations, ).
        num_iterations (int): number of iterations
                              (ie number of weak classifiers).

    Attributes:
        Xtrain (numpy.array): Array of flattened images (float32).
        ytrain (numpy.array): Labels array (float32).
        num_iterations (int): Number of iterations for the boosting loop.
        weakClassifiers (list): List of weak classifiers appended in each
                               iteration.
        alphas (list): List of alpha values, one for each classifier.
        num_obs (int): Number of observations.
        weights (numpy.array): Array of normalized weights, one for each
                               observation.
        eps (float): Error threshold value to indicate whether to update
                     the current weights or stop training.
    """

    def __init__(self, X, y, num_iterations):
        self.Xtrain = np.float32(X)
        self.ytrain = np.float32(y)
        self.num_iterations = num_iterations
        self.weakClassifiers = []
        self.alphas = []
        self.num_obs = X.shape[0]
        self.weights = np.array([1.0 / self.num_obs] * self.num_obs)  # uniform weights
        self.eps = 0.0001

    def train(self):
        """Implement the for loop shown in the problem set instructions."""
        for j in range(0, self.num_iterations):
            h = WeakClassifier(X=self.Xtrain, y=self.ytrain, weights=self.weights)
            h.train()
            mod_j = h.predict(np.transpose(self.Xtrain))
            inds = self.ytrain != mod_j
            ej = np.sum(self.weights[inds])/np.sum(self.weights)
            a = 0.5 * np.log((1. - ej)/ej)
            self.weakClassifiers = np.append(self.weakClassifiers, h)
            self.alphas = np.append(self.alphas, a)
            if ej > self.eps:
                self.weights[inds] = self.weights[inds] * np.exp(-a * mod_j[inds] * self.ytrain[inds])
            else:
                break


    def evaluate(self):
        """Return the number of correct and incorrect predictions.

        Use the training data (self.Xtrain) to obtain predictions. Compare
        them with the training labels (self.ytrain) and return how many
        where correct and incorrect.

        Returns:
            tuple: two-element tuple containing:
                correct (int): Number of correct predictions.
                incorrect (int): Number of incorrect predictions.
        """
        pred = self.predict(self.Xtrain)
        grdt  = self.ytrain
        correct = np.sum(grdt==pred)
        incorrect = np.sum(grdt!=pred)
        return (correct, incorrect)

    def predict(self, X):
        """Return predictions for a given array of observations.

        Use the alpha values stored in self.aphas and the weak classifiers
        stored in self.weakClassifiers.

        Args:
            X (numpy.array): Array of flattened images (observations).

        Returns:
            numpy.array: Predictions, one for each row in X.
        """
        H = [[m.predict(np.transpose(X))] for m in self.weakClassifiers]
        for i in range(0, len(self.alphas)):
            H[i] = self.alphas[i] * np.array(H[i])
        H = np.sum(H, axis=0)
        H = H[0]
        H = np.sign(H)
        return H


class HaarFeature:
    """Haar-like features.

    Args:
        feat_type (tuple): Feature type {(2, 1), (1, 2), (3, 1), (2, 2)}.
        position (tuple): (row, col) position of the feature's top left corner.
        size (tuple): Feature's (height, width)

    Attributes:
        feat_type (tuple): Feature type.
        position (tuple): Feature's top left corner.
        size (tuple): Feature's width and height.
    """

    def __init__(self, feat_type, position, size):
        self.feat_type = feat_type
        self.position = position
        self.size = size

    def _create_two_horizontal_feature(self, shape):
        """Create a feature of type (2, 1).

        Use int division to obtain half the height.

        Args:
            shape (tuple):  Array numpy-style shape (rows, cols).

        Returns:
            numpy.array: Image containing a Haar feature. (uint8).
        """
        f = np.zeros(shape)
        y, x = self.position
        h, w = self.size
        f[y:y+int(h/2), x:x+w] = 255
        f[y+int(h/2):y+h, x:x+w] = 126
        return f

    def _create_two_vertical_feature(self, shape):
        """Create a feature of type (1, 2).

        Use int division to obtain half the width.

        Args:
            shape (tuple):  Array numpy-style shape (rows, cols).

        Returns:
            numpy.array: Image containing a Haar feature. (uint8).
        """
        f = np.zeros(shape)
        y, x = self.position
        h, w = self.size
        f[y:y+h, x:x+int(w/2)] = 255
        f[y:y+h, x+int(w/2):x+w] = 126
        return f

    def _create_three_horizontal_feature(self, shape):
        """Create a feature of type (3, 1).

        Use int division to obtain a third of the height.

        Args:
            shape (tuple):  Array numpy-style shape (rows, cols).

        Returns:
            numpy.array: Image containing a Haar feature. (uint8).
        """
        f = np.zeros(shape)
        y, x = self.position
        h, w = self.size
        f[y:y+int(h/3), x:x+w] = 255
        f[y+int(h/3):y+int(h/3)*2, x:x+w] = 126
        f[y+int(h/3)*2:y+h, x:x+w] = 255
        return f

    def _create_three_vertical_feature(self, shape):
        """Create a feature of type (1, 3).

        Use int division to obtain a third of the width.

        Args:
            shape (tuple):  Array numpy-style shape (rows, cols).

        Returns:
            numpy.array: Image containing a Haar feature. (uint8).
        """
        f = np.zeros(shape)
        y, x = self.position
        h, w = self.size
        f[y:y+h, x:x+int(w/3)] = 255
        f[y:y+h, x+int(w/3):x+int(2*w/3)] = 126
        f[y:y+h, x+int(2*w/3):x+w] = 255
        return f

    def _create_four_square_feature(self, shape):
        """Create a feature of type (2, 2).

        Use int division to obtain half the width and half the height.

        Args:
            shape (tuple):  Array numpy-style shape (rows, cols).

        Returns:
            numpy.array: Image containing a Haar feature. (uint8).
        """
        f = np.zeros(shape)
        y, x = self.position
        h, w = self.size
        f[y:y+int(h/2), x:x+int(w/2)] = 126
        f[y:y+int(h/2), x+int(w/2):x+w] = 255
        f[y+int(h/2):y+h, x:x+int(w/2)] = 255
        f[y+int(h/2):y+h, x+int(w/2):x+w] = 126
        return f

    def preview(self, shape=(24, 24), filename=None):
        """Return an image with a Haar-like feature of a given type.

        Function that calls feature drawing methods. Each method should
        create an 2D zeros array. Each feature is made of a white area (255)
        and a gray area (126).

        The drawing methods use the class attributes position and size.
        Keep in mind these are in (row, col) and (height, width) format.

        Args:
            shape (tuple): Array numpy-style shape (rows, cols).
                           Defaults to (24, 24).

        Returns:
            numpy.array: Array containing a Haar feature (float or uint8).
        """

        if self.feat_type == (2, 1):  # two_horizontal
            X = self._create_two_horizontal_feature(shape)

        if self.feat_type == (1, 2):  # two_vertical
            X = self._create_two_vertical_feature(shape)

        if self.feat_type == (3, 1):  # three_horizontal
            X = self._create_three_horizontal_feature(shape)

        if self.feat_type == (1, 3):  # three_vertical
            X = self._create_three_vertical_feature(shape)

        if self.feat_type == (2, 2):  # four_square
            X = self._create_four_square_feature(shape)

        if filename is None:
            cv2.imwrite("output/{}_feature.png".format(self.feat_type), X)

        else:
            cv2.imwrite("output/{}.png".format(filename), X)

        return X

    def evaluate(self, ii):
        """Evaluates a feature's s on a given integral image.

        Calculate the s of a feature defined by the self.feat_type.
        Using the integral image and the sum / subtraction of rectangles to
        obtain a feature's value. Add the feature's white area value and
        subtract the gray area.

        For example, on a feature of type (2, 1):
        s = sum of pixels in the white area - sum of pixels in the gray area

        Keep in mind you will need to use the rectangle sum / subtraction
        method and not numpy.sum(). This will make this process faster and
        will be useful in the ViolaJones algorithm.

        Args:
            ii (numpy.array): Integral Image.

        Returns:
            float: Score value.
        """
        h, w = self.size
        y, x = self.position[:2]
        y1, x1 = self.position[:2]

        if self.feat_type == (2, 1):
            y2, x2 = y+int(h/2)-1,  x1+w-1
            a = ii[y2, x2]-ii[y1-1, x2]-ii[y2, x1-1]+ii[y1-1, x1-1]
            y1 = y+int(h/2)
            x1 = x
            y2, x2 = y+int(h)-1, x1+w-1
            b = ii[y2, x2]-ii[y1-1, x2]-ii[y2, x1-1]+ii[y1-1, x1-1]
            s = a-b

        if self.feat_type == (1, 2):
            y2, x2 = y1+int(h)-1,  x1+int(w/2)-1
            a = ii[y2, x2]-ii[y1-1, x2]-ii[y2, x1-1]+ii[y1-1, x1-1]
            y1 = y
            x1 = x+int(w/2)
            y2, x2 =y1+int(h)-1, x+int(w)-1
            b = ii[y2, x2]-ii[y1-1, x2]-ii[y2, x1-1]+ii[y1-1, x1-1]
            s = a-b

        if self.feat_type == (3, 1):
            y2 = y+int(h/3)-1
            x2 = x1+w-1
            a = ii[y2, x2]-ii[y1-1, x2]-ii[y2, x1-1]+ii[y1-1, x1-1]
            y1 = y+int(h/3)
            x1 = x
            y2, x2 = y+int(h/3)+int(h/3)-1, x1+w-1
            b = ii[y2, x2]-ii[y1-1, x2]-ii[y2, x1-1]+ii[y1-1, x1-1]
            y1 = y+int(h/3)+int(h/3)
            x1 = x
            y2, x2 = y+h-1, x1+w-1
            c = ii[y2, x2]-ii[y1-1,x2]-ii[y2,x1-1]+ii[y1-1,x1-1]
            s = a-b+c

        if self.feat_type == (1, 3):
            y2 = y1+int(h)-1
            x2 = x1+int(w/3)-1
            a = ii[y2, x2]-ii[y1-1, x2]-ii[y2, x1-1]+ii[y1-1, x1-1]
            y1 = y
            x1 = x+int(w/3)
            y2, x2 = y1+int(h)-1, x+int(w/3)+int(w/3)-1
            b = ii[y2, x2]-ii[y1-1, x2]-ii[y2, x1-1]+ii[y1-1, x1-1]
            y1 = y
            x1 = x+int(w/3)+int(w/3)
            y2, x2 = y+h-1, x+int(w)-1
            c = ii[y2, x2]-ii[y1-1, x2]-ii[y2, x1-1]+ii[y1-1, x1-1]
            s = a-b+c

        if self.feat_type == (2, 2):
            y2 = y1+int(h/2)-1
            x2 = x1+int(w/2)-1
            a = ii[y2, x2]-ii[y1-1,x2]-ii[y2,x1-1]+ii[y1-1,x1-1]
            y1 = y
            x1 = x+int(w/2)
            y2, x2 = y1+int(h/2)-1, x+w-1
            b = ii[y2, x2]-ii[y1-1, x2]-ii[y2, x1-1]+ii[y1-1, x1-1]
            y1 = y+int(h/2)
            x1 = x
            y2, x2 = y+h-1, x1+int(w/2)-1
            c = ii[y2, x2]-ii[y1-1, x2]-ii[y2, x1-1]+ii[y1-1, x1-1]
            y1 = y+int(h/2)
            x1 = x+int(w/2)
            y2, x2 = y+h-1, x+w-1
            d = ii[y2, x2]-ii[y1-1, x2]-ii[y2, x1-1]+ii[y1-1, x1-1]
            s = -d-a+b+c

        return s


def convert_images_to_integral_images(images):
    """Convert a list of grayscale images to integral images.

    Args:
        images (list): List of grayscale images (uint8 or float).

    Returns:
        (list): List of integral images.
    """
    lst = []
    for i in images:
        lst.append(np.cumsum(np.cumsum(i, axis=0), axis=1))
    return lst


class ViolaJones:
    """Viola Jones face detection method

    Args:
        pos (list): List of positive images.
        neg (list): List of negative images.
        integral_images (list): List of integral images.

    Attributes:
        haarFeatures (list): List of haarFeature objects.
        integralImages (list): List of integral images.
        classifiers (list): List of weak classifiers (VJ_Classifier).
        alphas (list): Alpha values, one for each weak classifier.
        posImages (list): List of positive images.
        negImages (list): List of negative images.
        labels (numpy.array): Positive and negative labels.
    """
    def __init__(self, pos, neg, integral_images):
        self.haarFeatures = []
        self.integralImages = integral_images
        self.classifiers = []
        self.alphas = []
        self.posImages = pos
        self.negImages = neg
        self.labels = np.hstack((np.ones(len(pos)), -1*np.ones(len(neg))))
        self.threshold = 1.0

    def createHaarFeatures(self):
        # Let's take detector resolution of 24x24 like in the paper
        FeatureTypes = {"two_horizontal": (2, 1),
                        "two_vertical": (1, 2),
                        "three_horizontal": (3, 1),
                        "three_vertical": (1, 3),
                        "four_square": (2, 2)}

        haarFeatures = []
        for _, feat_type in FeatureTypes.items():
            for sizei in range(feat_type[0], 24 + 1, feat_type[0]):
                for sizej in range(feat_type[1], 24 + 1, feat_type[1]):
                    for posi in range(0, 24 - sizei + 1, 4):
                        for posj in range(0, 24 - sizej + 1, 4):
                            haarFeatures.append(
                                HaarFeature(feat_type, [posi, posj],
                                            [sizei-1, sizej-1]))
        self.haarFeatures = haarFeatures

    def set_threshold(self, threshold):
        self.threshold = threshold

    def init_train(self):
        """ This function initializes self.scores, self.weights

        Args:
            None

        Returns:
            None
        """
    
        # Use this scores array to train a weak classifier using VJ_Classifier
        # in the for loop below.
        if not self.integralImages or not self.haarFeatures:
            print("No images provided. run convertImagesToIntegralImages() first")
            print("       Or no features provided. run creatHaarFeatures() first")
            return

        self.scores = np.zeros((len(self.integralImages), len(self.haarFeatures)))
        print(" -- compute all scores --")
        for i, im in enumerate(self.integralImages):
            self.scores[i, :] = [hf.evaluate(im) for hf in self.haarFeatures]

        weights_pos = np.ones(len(self.posImages), dtype='float') * 1.0 / (
                           2*len(self.posImages))
        weights_neg = np.ones(len(self.negImages), dtype='float') * 1.0 / (
                           2*len(self.negImages))
        self.weights = np.hstack((weights_pos, weights_neg))

    def train(self, num_classifiers):
        """ Initialize and train Viola Jones face detector

        The function should modify self.weights, self.classifiers, self.alphas, and self.threshold

        Args:
            None

        Returns:
            None
        """
        self.init_train()
        print(" -- select classifiers --")
        for t in range(num_classifiers):
            self.weights = self.weights/np.sum(self.weights)
            hj = VJ_Classifier(self.scores, self.labels, self.weights)
            hj.train()
            self.classifiers.append(hj)
            b = hj.error/(1-hj.error)
            a = np.log(1/b)
            for i in range(0, len(self.weights)):
                if hj.predict(self.scores[i]) == self.labels[i]:
                    ei = 0
                else:
                    ei = 1
                self.weights[i] =self.weights[i]*(b**(1-ei))
            self.alphas.append(a)

    def predict(self, images):
        """Return predictions for a given list of images.

        Args:
            images (list of element of type numpy.array): list of images (observations).

        Returns:
            list: Predictions, one for each element in images.
        """

        ii = convert_images_to_integral_images(images)

        scores = np.zeros((len(ii), len(self.haarFeatures)))

        # Populate the score location for each classifier 'clf' in
        # self.classifiers.
        for clf in self.classifiers:
        # Obtain the Haar feature id from clf.feature
            fid = clf.feature

        # Use this id to select the respective feature object from
        # self.haarFeatures
            hf = self.haarFeatures[fid]

        # Add the score value to score[x, feature id] calling the feature's
        # evaluate function. 'x' is each image in 'ii'

            for x, img in enumerate(ii):
                scores[x, fid] = hf.evaluate(img)

        result = []

        # Append the results for each row in 'scores'. This value is obtained
        # using the equation for the strong classifier H(x).

        for x in scores:
            H = 0
            t = 0
            for clf in self.classifiers:
                H = H+clf.predict(x)*self.alphas[t]
                t = t+1
            if H >= 0.5*np.sum(np.array(self.alphas)):
                result.append(1)
            else:
                result.append(-1)
        return result

    def faceDetection(self, image, filename):
        """Scans for faces in a given image.

        Complete this function following the instructions in the problem set
        document.

        Use this function to also save the output image.

        Args:
            image (numpy.array): Input image.
            filename (str): Output image file name.

        Returns:
            None.
        """
        frame = image.copy()
        img = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
        w = 24
        ptch = []
        ul = []
        lr = []
        for x in range(0,img.shape[0]-w):
            for y in range(0,img.shape[1]-w):
                ul.append([y, x])
                lr.append([y+w, x+w])
                ptch.append(np.array(img[x:x+w, y:y+w]))
        pred = self.predict(ptch)
        
        ulm = (np.array(ul)[np.array(pred) == 1]).mean(axis=0)
        lrm = (np.array(lr)[np.array(pred) == 1]).mean(axis=0)
        cv2.rectangle(frame, tuple(ulm.astype(int)), tuple(lrm.astype(int)), (255, 0, 0), 1)
        cv2.imwrite(os.path.join("output", filename + '.jpg'), frame)


class CascadeClassifier:
    """Viola Jones Cascade Classifier Face Detection Method

    Lesson: 8C-L2, Boosting and face detection

    Args:
        f_max (float): maximum acceptable false positive rate per layer
        d_min (float): minimum acceptable detection rate per layer
        f_target (float): overall target false positive rate
        pos (list): List of positive images.
        neg (list): List of negative images.

    Attributes:
        f_target: overall false positive rate
        classifiers (list): Adaboost classifiers
        train_pos (list of numpy arrays):  
        train_neg (list of numpy arrays): 

    """
    def __init__(self, pos, neg, f_max_rate=0.30, d_min_rate=0.70, f_target = 0.07):
        
        train_percentage = 0.85

        pos_indices = np.random.permutation(len(pos)).tolist()
        neg_indices = np.random.permutation(len(neg)).tolist()

        train_pos_num = int(train_percentage * len(pos))
        train_neg_num = int(train_percentage * len(neg))

        pos_train_indices = pos_indices[:train_pos_num]
        pos_validate_indices = pos_indices[train_pos_num:]

        neg_train_indices = neg_indices[:train_neg_num]
        neg_validate_indices = neg_indices[train_neg_num:]

        self.train_pos = [pos[i] for i in pos_train_indices]
        self.train_neg = [neg[i] for i in neg_train_indices]

        self.validate_pos = [pos[i] for i in pos_validate_indices]
        self.validate_neg = [neg[i] for i in neg_validate_indices]

        self.f_max_rate = f_max_rate
        self.d_min_rate = d_min_rate
        self.f_target = f_target
        self.classifiers = []

    def predict(self, classifiers, img):
        """Predict face in a single image given a list of cascaded classifiers

        Args:
            classifiers (list of element type ViolaJones): list of ViolaJones classifiers to predict 
                where index i is the i'th consecutive ViolaJones classifier
            img (numpy.array): Input image

        Returns:
            Return 1 (face detected) or -1 (no face detected) 
        """

        # TODO
        raise NotImplementedError

    def evaluate_classifiers(self, pos, neg, classifiers):
        """ 
        Given a set of classifiers and positive and negative set
        return false positive rate and detection rate 

        Args:
            pos (list): Input image.
            neg (list): Output image file name.
            classifiers (list):  

        Returns:
            f (float): false positive rate
            d (float): detection rate
            false_positives (list): list of false positive images
        """

        # TODO
        raise NotImplementedError

    def train(self):
        """ 
        Trains a cascaded face detector

        Sets self.classifiers (list): List of ViolaJones classifiers where index i is the i'th consecutive ViolaJones classifier

        Args:
            None

        Returns:
            None
             
        """
        # TODO
        raise NotImplementedError


    def faceDetection(self, image, filename="ps6-5-b-1.jpg"):
        """Scans for faces in a given image using the Cascaded Classifier.

        Complete this function following the instructions in the problem set
        document.

        Use this function to also save the output image.

        Args:
            image (numpy.array): Input image.
            filename (str): Output image file name.

        Returns:
            None.
        """
        raise NotImplementedError
