<img src="hisat_dist_plot_before_diffexp.png" alt="hisat distribution plot before normalization" style="height: 500px; width:500px;"/><br>
This is a distribution plot of our samples using the hisat alignment method before noramlization. There are not any clearly visable patterns.<br>
<br>
<br>
<img src="hisat_pca_before_deseq2.png" alt="hisat pca plot before normalization" style="height: 500px; width:500px;"/><br>
This is a pca plot of our samples using the hisat alignment method before normalization. There is clear separation between the samples from the liver and intestines. This means the pca after normalization should look good as well.<br>
<br>
<br>
<img src="hisat_pca_after_deseq2.png" alt="hisat pca plot after normalization" style="height: 500px; width:500px;"/><br>
This is a pca plot of our samples using the hisat alignment method after normalization with DESeq2. There is clear clustering of the samples from the liver and from the intestines. This leads me to believe that there is a signifficant different between miracidia from the liver and intestines.<br>
<br>
<br>
<img src="hisat_volcano_plot.png" alt="hisat volcano plot after normalization" style="height: 500px; width:500px;"/><br>
This is a volcano plot of our samples using the hisat alignment method after normalization with DESeq2. Most genes had a padj value of 0. This could be due to using from the reduced dataset rather than the whole dataset as there is less variation.<br>
<br>
<br>
<img src="star_dist_plot_before_diffexp.png" alt="star distribution plot before normalization" style="height: 500px; width:500px;"/><br>
This is a distribution plot of our samples using the star alignment method before noramlization. There are not any clearly visable patterns.<br>
<br>
<br>
<img src="star_pca_before_deseq2.png" alt="star pca plot before normalization" style="height: 500px; width:500px;"/><br>
This is a pca plot of our samples using the star alignment method before normalization. There is clear separation between the samples from the liver and intestines. This means the pca after normalization should look good as well.<br>
<br>
<br>
<img src="star_pca_after_deseq2.png" alt="star pca plot after normalization" style="height: 500px; width:500px;"/><br>
This is a pca plot of our samples using the star alignment method after normalization with DESeq2. There is clear clustering of the samples from the liver and from the intestines. This leads me to believe that there is a signifficant different between miracidia from the liver and intestines.<br>
<br>
<br>
<img src="star_volcano_plot.png" alt="star volcano plot after normalization" style="height: 500px; width:500px;"/><br>
This is a volcano plot of our samples using the star alignment method after normalization with DESeq2. Most genes had a padj value of 0. This could be due to using from the reduced dataset rather than the whole dataset as there is less variation.<br>
<br>
<br>
Overall, the plots show that the miracidia from the liver are transcriptomically different from those in the intestines.