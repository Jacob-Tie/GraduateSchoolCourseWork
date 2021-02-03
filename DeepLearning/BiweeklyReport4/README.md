# Biweekly Progress Report 4
# Jacob Tiede
## File Descriptions:
### ImageToImageTranslation.ipynb
In this file I will summarize and implement a paper that seeks to perform image to image translation on unlabeled data. This is specifically applied to taking an image taken at a particular time of day and translating it into the same picture taken at a different time of the day.
### LSTMsRevisited.ipynb
Here I address feedback from my prior biweekly report related to LSTMs. I speed up the training process by increasing the batch size, though this seems to encourage overfitting. I also look at some of the mistakes that my LSTM makes and I implement a spell checker, which did not have a substantial effect on accuracy.
### ProjectionGANOnMNIST.ipynb
In this file I use one of the concepts that I detailed in ImageToImageTranslation.ipynb (projection GANs) on the MNIST dataset. I compare using a projection GAN vs using a more conventional conditional GAN, neither of which seemed to converge very quickly. To try and speed up the training process I implemented an idea of my own where I use a pretrained image classification neural net to determine what digit is being generated, which I used as part of the loss for the neural net. This substantially improved result in creating images of a given class
## Final Thoughts on this Progress Report
I think that I've learned many useful skills this week. This included how summary statistics can be used to glean information about images, as well as different perspectives on how a good GAN should be structured. I was also able to successfully implement one of my own ideas (perhaps it's already in the literature, but if it is I haven't found it) which actually seemed to help the training process for conditional GANs.
