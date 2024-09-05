
# Facilitating precision medicine in HCC patients by DL-directed lncRNAs classification and ascertaining causal markers

Project Description:
The project created a Deep Learning (DL) based hierarchical composite framework comprising three Deep Neural Network (DNN) models. It is able to classify hepatocellular carcinoma samples from normal to the tumor, and further as Early (Stage I) and Advanced, also distinguished as Stage II, Stage III and Stage IV, based on RNA-seq expression values of long non-coding RNAs (lncRNAs). This aided in biomarker discovery.  

User Manual: 
This user manual contains the steps required to load and run the project in Python. 

1. Installing Spyder 
Any Python IDE can be used. Here the operation with Spyder is explained. Spyder IDE can be used online with https://mybinder.org/, installed using external distributor, e.g. Anaconda or stand-alone in your system. Refer to https://docs.spyder-ide.org/current/installation.html to install the latest version of Spyder as per your system. Create a new environment in conda as depicted in the HTML page.

2. It is always recommended to set directory and download all files there for ease of use. This can be done by:  
pip install os  
import os  
os.chdir(path)   

  In the last command, in place of path, you can add the path of your directory. 

2. After this, install the necessary libraries to run the model. These include rpy2, numpy, pandas, scipy.stats, warnings, keras and tkinter. The commands for this are:   
pip install rpy2  
pip install numpy  
pip install pandas  
pip install scipy.stats  
pip install warnings  
pip install keras  
pip install tkinter  

  As you type these commands one by one, press Enter and it will start installing. It might take some time but once it is done, you are ready to run the model. 

3. Download the data file in CSV (comma-separated) format to submit as input when asked while running the model. 

4. Open the script file in Spyder and click on run file to run the program. You will be asked to upload the data once the code starts running. Upload the downloaded data file and let the code run to completion.  

5. An output file Results.csv will automatically be downloaded in your directory. It can be opened in MS Excel or Google Sheets and the classification results can be interpreted for each Patient ID either as Normal or Tumor, Early or Advanced Tumor or AJCC Pathological Stage. 

