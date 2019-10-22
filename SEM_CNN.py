from keras.models import Sequential
from keras.preprocessing.image import ImageDataGenerator, img_to_array, load_img
from keras.layers import Dense, Activation, Flatten, Dropout, BatchNormalization
from keras.layers import Conv2D, MaxPooling2D
from keras import regularizers, optimizers
import pandas as pd
import numpy as np
from PIL import Image
import os
from sklearn.utils import shuffle
from keras.utils import to_categorical


df = pd.read_csv("/Users/nbense/NASA/BRAILLE/SEM_CNN/SEM_CNN_Labels.csv")
df["labels"]=df["labels"].apply(lambda x:x.split(","))
df = shuffle(df)

datagen=ImageDataGenerator(rescale=1./255.)
test_datagen=ImageDataGenerator(rescale=1./255.)

train_generator=datagen.flow_from_dataframe(
dataframe=df[:600],
directory="/Users/nbense/NASA/BRAILLE/SEM_CNN/Images/",
x_col="Filenames",
y_col="labels",
batch_size=32,
seed=42,
shuffle=True,
class_mode="categorical",
classes=["filament", "film", "spheroid", "fuzzy", "overview", "smooth", "none", "rod", "mineral"],
target_size=(1020,1024))

valid_generator=test_datagen.flow_from_dataframe(
dataframe=df[600:670],
directory="/Users/nbense/NASA/BRAILLE/SEM_CNN/Images/",
x_col="Filenames",
y_col="labels",
batch_size=32,
seed=42,
shuffle=True,
class_mode="categorical",
classes=["filament", "film", "spheroid", "fuzzy", "overview", "smooth", "none", "rod", "mineral"],
target_size=(1020,1024))

test_generator=test_datagen.flow_from_dataframe(
dataframe=df[670:733],
directory="/Users/nbense/NASA/BRAILLE/SEM_CNN/Images/",
x_col="Filenames",
batch_size=1,
seed=42,
shuffle=False,
class_mode=None,
target_size=(1020,1024))

model = Sequential()
model.add(Conv2D(32, (3, 3), padding='same',
                 input_shape=(1020,1024,3)))
model.add(Activation('relu'))
model.add(Conv2D(32, (3, 3)))
model.add(Activation('relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))
model.add(Dropout(0.25))
model.add(Conv2D(64, (3, 3), padding='same'))
model.add(Activation('relu'))
model.add(Conv2D(64, (3, 3)))
model.add(Activation('relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))
model.add(Dropout(0.25))
model.add(Flatten())
model.add(Dense(512))
model.add(Activation('relu'))
model.add(Dropout(0.5))
model.add(Dense(9, activation='sigmoid'))
model.compile(optimizers.rmsprop(lr=0.0001, decay=1e-6),loss="binary_crossentropy",metrics=["accuracy"])

STEP_SIZE_TRAIN=train_generator.n//train_generator.batch_size
STEP_SIZE_VALID=valid_generator.n//valid_generator.batch_size
STEP_SIZE_TEST=test_generator.n//test_generator.batch_size
model.fit_generator(generator=train_generator,
                    steps_per_epoch=STEP_SIZE_TRAIN,
                    validation_data=valid_generator,
                    validation_steps=STEP_SIZE_VALID,
                    epochs=10
)

test_generator.reset()
pred=model.predict_generator(test_generator,
steps=STEP_SIZE_TEST,
verbose=1)

pred_bool = (pred >0.5)
predictions=[]
labels = train_generator.class_indices
labels = dict((v,k) for k,v in labels.items())
for row in pred_bool:
    l=[]
    for index,cls in enumerate(row):
        if cls:
            l.append(labels[index])
    predictions.append(",".join(l))
filenames=test_generator.filenames
results=pd.DataFrame({"Filename":filenames,
                      "Predictions":predictions})

results.to_csv("results.csv",index=False)
