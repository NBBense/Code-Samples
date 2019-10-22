from sklearn.decomposition import PCA
import pandas as pd
from PIL import Image
import numpy as np
import os
import matplotlib
from sklearn.manifold import TSNE
matplotlib.use('TkAgg')

filenames = []
for file in os.listdir('/Users/nbense/NASA/BRAILLE/SEM_Updated_All/tif_files/'):
    if file.endswith( ('.tif') ): 
        print(file)
        filenames.append(file)
with open("/Users/nbense/NASA/BRAILLE/SE_Subset.txt", "r") as subsetfile:
    for newline in subsetfile:
        filenames.append(newline.rstrip())

os.chdir('/Users/nbense/NASA/BRAILLE/SEM_Updated_All/tif_files/')
img1 = Image.open(str(filenames[0])).convert("L")
listsamples = []
listfilenames = []
listsamples.append(filenames[0].split('-')[0])
listfilenames.append(filenames[0])
print (listsamples)
imgarr1=np.array(img1)
imgarr1=imgarr1[0:1020].flatten()
img2 = Image.open(str(filenames[1])).convert("L")
listsamples.append(filenames[1].split('-')[0])
listfilenames.append(filenames[1])
imgarr2=np.array(img2)
imgarr2=imgarr2[0:1020].flatten()
finalimgarr=np.append([imgarr1],[imgarr2], axis=0)
for filename in filenames[2:]:
    fields = filename.split('-')
    print ((fields[0])[:3])
    img = Image.open(str(filename)).convert("L")
    imgarr = np.array(img)
    imgarr = imgarr[0:1020].flatten()
    print (filename)
    listsamples.append((fields[0])[:3])
    listfilenames.append(filename)
    print len(imgarr)
    print len(finalimgarr)
    finalimgarr = np.append(finalimgarr,[imgarr], axis=0)

feat_cols = [ 'pixel'+str(i) for i in range(finalimgarr.shape[1]) ]
df = pd.DataFrame(finalimgarr,columns=feat_cols)
imglabel = np.array(listsamples, dtype=object)
df['label'] = imglabel
namefile = np.array(listfilenames, dtype=object)
df['filename'] = namefile

pca = PCA(n_components=3)
pca_result = pca.fit_transform(df[feat_cols].values)

df['pca-one'] = pca_result[:,0]
df['pca-two'] = pca_result[:,1] 
df['pca-three'] = pca_result[:,2]

print ('Explained variation per principal component: {}'.format(pca.explained_variance_ratio_))

pca_50 = PCA(n_components=50)
pca_result_50 = pca_50.fit_transform(df[feat_cols].values)
print ('Cumulative explained variation for 50 principal components: {}'.format(np.sum(pca_50.explained_variance_ratio_)))

new_pca_df = df.loc[:, ['pca-one','pca-two','pca-three','label','filename']]
new_pca_df.to_csv('/Users/nbense/NASA/BRAILLE/SEM_Updated_All/SEM_PCA_SE_Subset_Labels.csv')

tsne = TSNE(n_components=2, verbose=1, perplexity=40, n_iter=300)
tsne.fit_transform(df[feat_cols].values)

tsne_df = df[feat_cols].copy()
tsne_results = tsne.fit_transform(df[feat_cols].values)
tsne_df['x-tsne'] = tsne_results[:,0]
tsne_df['y-tsne'] = tsne_results[:,1]
tsne_df['label'] = imglabel
tsne_df['filename'] = namefile

new_tsne_df = tsne_df.loc[:, ['x-tsne','y-tsne','label','filename']]
new_tsne_df.to_csv('/Users/nbense/NASA/BRAILLE/SEM_Updated_All/SEM_TSNE_SE_Subset_Labels.csv')
