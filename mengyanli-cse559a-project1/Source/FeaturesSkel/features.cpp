#include <assert.h>
#include <math.h>
#include <FL/Fl.H>
#include <FL/Fl_Image.H>
#include "features.h"
#include "ImageLib/FileIO.h"

#define PI 3.14159265358979323846

// Compute features of an image.
bool computeFeatures(CFloatImage &image, FeatureSet &features, int featureType) {
	// TODO: Instead of calling dummyComputeFeatures, write your own
	// feature computation routines and call them here.
	switch (featureType) {
	case 1:
		dummyComputeFeatures(image, features);
		break;
	case 2:
		ComputeHarrisFeatures(image, features);
		break;
	default:
		return false;
	}

	// This is just to make sure the IDs are assigned in order, because
	// the ID gets used to index into the feature array.
	for (unsigned int i=0; i<features.size(); i++) {
		features[i].id = i+1;
	}

	return true;
}

// Perform a query on the database.  This simply runs matchFeatures on
// each image in the database, and returns the feature set of the best
// matching image.
bool performQuery(const FeatureSet &f, const ImageDatabase &db, int &bestIndex, vector<FeatureMatch> &bestMatches, double &bestScore, int matchType) {
	// Here's a nice low number.
	bestScore = -1e100;

	vector<FeatureMatch> tempMatches;
	double tempScore;

	for (unsigned int i=0; i<db.size(); i++) {
		if (!matchFeatures(f, db[i].features, tempMatches, tempScore, matchType)) {
			return false;
		}

		if (tempScore > bestScore) {
			bestIndex = i;
			bestScore = tempScore;
			bestMatches = tempMatches;
		}
	}

	return true;
}

// Match one feature set with another.
bool matchFeatures(const FeatureSet &f1, const FeatureSet &f2, vector<FeatureMatch> &matches, double &totalScore, int matchType) {
	// TODO: We have given you the ssd matching function, you must write your own
	// feature matching function for the ratio test.
	
	printf("\nMatching features.......\n");

	switch (matchType) {
	case 1:
		ssdMatchFeatures(f1, f2, matches, totalScore);
		return true;
	case 2:
		ratioMatchFeatures(f1, f2, matches, totalScore);
		return true;
	default:
		return false;
	}
}

// Evaluate a match using a ground truth homography.  This computes the
// average SSD distance between the matched feature points and
// the actual transformed positions.
double evaluateMatch(const FeatureSet &f1, const FeatureSet &f2, const vector<FeatureMatch> &matches, double h[9]) {
	double d = 0;
	int n = 0;

	double xNew;
	double yNew;

    unsigned int num_matches = matches.size();
	for (unsigned int i=0; i<num_matches; i++) {
		int id1 = matches[i].id1;
        int id2 = matches[i].id2;
        applyHomography(f1[id1-1].x, f1[id1-1].y, xNew, yNew, h);
		d += sqrt(pow(xNew-f2[id2-1].x,2)+pow(yNew-f2[id2-1].y,2));
		n++;
	}	

	return d / n;
}

void addRocData(const FeatureSet &f1, const FeatureSet &f2, const vector<FeatureMatch> &matches, double h[9],vector<bool> &isMatch,double threshold,double &maxD) {
	double d = 0;

	double xNew;
	double yNew;

    unsigned int num_matches = matches.size();
	for (unsigned int i=0; i<num_matches; i++) {
		int id1 = matches[i].id1;
        int id2 = matches[i].id2;
		applyHomography(f1[id1-1].x, f1[id1-1].y, xNew, yNew, h);

		// Ignore unmatched points.  There might be a better way to
		// handle this.
		d = sqrt(pow(xNew-f2[id2-1].x,2)+pow(yNew-f2[id2-1].y,2));
		if (d<=threshold)
		{
			isMatch.push_back(1);
		}
		else
		{
			isMatch.push_back(0);
		}

		if (matches[i].score>maxD)
			maxD=matches[i].score;
	}	
}

vector<ROCPoint> computeRocCurve(vector<FeatureMatch> &matches,vector<bool> &isMatch,vector<double> &thresholds)
{
	vector<ROCPoint> dataPoints;

	for (int i=0; i < (int)thresholds.size();i++)
	{
		//printf("Checking threshold: %lf.\r\n",thresholds[i]);
		int tp=0;
		int actualCorrect=0;
		int fp=0;
		int actualError=0;
		int total=0;

        int num_matches = (int) matches.size();
		for (int j=0;j < num_matches;j++)
		{
			if (isMatch[j])
			{
				actualCorrect++;
				if (matches[j].score<thresholds[i])
				{
					tp++;
				}
			}
			else
			{
				actualError++;
				if (matches[j].score<thresholds[i])
				{
					fp++;
				}
            }
			
			total++;
		}

		ROCPoint newPoint;
		//printf("newPoints: %lf,%lf",newPoint.trueRate,newPoint.falseRate);
		newPoint.trueRate=(double(tp)/actualCorrect);
		newPoint.falseRate=(double(fp)/actualError);
		//printf("newPoints: %lf,%lf",newPoint.trueRate,newPoint.falseRate);

		dataPoints.push_back(newPoint);
	}

	return dataPoints;
}


// Compute silly example features.  This doesn't do anything
// meaningful.
void dummyComputeFeatures(CFloatImage &image, FeatureSet &features) {
	CShape sh = image.Shape();
	Feature f;

	for (int y=0; y<sh.height; y++) {
		for (int x=0; x<sh.width; x++) {
			double r = image.Pixel(x,y,0);
			double g = image.Pixel(x,y,1);
			double b = image.Pixel(x,y,2);

			if ((int)(255*(r+g+b)+0.5) % 100  == 1) {
				// If the pixel satisfies this meaningless criterion,
				// make it a feature.
				
				f.type = 1;
				f.id += 1;
				f.x = x;
				f.y = y;

				f.data.resize(1);
				f.data[0] = r + g + b;

				features.push_back(f);
			}
		}
	}
}

void ComputeHarrisFeatures(CFloatImage &image, FeatureSet &features)
{
	//Create grayscale image used for Harris detection
	CFloatImage grayImage=ConvertToGray(image);

	//Create image to store Harris values
	CFloatImage harrisImage(image.Shape().width,image.Shape().height,1);

	//Create image to store local maximum harris values as 1, other pixels 0
	CByteImage harrisMaxImage(image.Shape().width,image.Shape().height,1);

	
	// CFloatImage orientationImage(image.Shape().width, image.Shape().height, 1);
	// getOrientation(grayImage);



	//compute Harris values puts harris values at each pixel position in harrisImage. 
	//You'll need to implement this function.
    computeHarrisValues(grayImage, harrisImage);
	
	// Threshold the harris image and compute local maxima.  
	//You'll need to implement this function.
	computeLocalMaxima(harrisImage,harrisMaxImage);


    
    // Prints out the harris image for debugging purposes
	CByteImage tmp(harrisImage.Shape());
	convertToByteImage(harrisImage, tmp);
    WriteFile(tmp, "harris.tga");
    

	// TO DO--------------------------------------------------------------------
	//Loop through feature points in harrisMaxImage and create feature descriptor 
	//for each point above a threshold

    int id = 0;
    for (int y = 0; y < harrisMaxImage.Shape().height; y++) {
		for (int x = 0; x < harrisMaxImage.Shape().width; x++) {
		
			// Skip over non-maxima
            if (harrisMaxImage.Pixel(x, y, 0) == 0)
                continue;

            //TO DO---------------------------------------------------------------------
		    // Fill in feature with descriptor data here. 
            Feature f;

            f.type = 2;
            f.id = id++;
            f.x = x;
            f.y = y;
            f.angleRadians = 0;

            // Add the feature to the list of features
            features.push_back(f);

         
        }
	}   
	//computeSimpleDescriptors(grayImage, features);
    computeMOPsDescriptors(grayImage, features);
}


//TO DO---------------------------------------------------------------------
//Loop through the image to compute the harris corner values as described in class
// srcImage:  grayscale of original image
// harrisImage:  populate the harris values per pixel in this image
void computeHarrisValues(CFloatImage &srcImage, CFloatImage &harrisImage)
{
    
	int w = srcImage.Shape().width;
    int h = srcImage.Shape().height;

    CFloatImage gaussianMask(5,5,1);    
    CFloatImage xGradient(w,h,1);    
    CFloatImage yGradient(w,h,1);    

 
    // generate 5x5 mask
    for (int i = 0; i < 5; i++)
    {
        for(int j = 0; j < 5; j++)
        {
            gaussianMask.Pixel(i,j,0) = gaussian5x5[i * 5 + j];
        }
    }
 
      
    //compute the gradient
    Convolve(srcImage,xGradient,ConvolveKernel_SobelX);
	Convolve(srcImage,yGradient,ConvolveKernel_SobelY);


	CFloatImage H_a(srcImage.Shape());
	CFloatImage H_b(srcImage.Shape());
	CFloatImage H_c(srcImage.Shape());


    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {  
        	H_a.Pixel(x,y,0) = xGradient.Pixel(x,y,0) * xGradient.Pixel(x,y,0);
			H_b.Pixel(x,y,0) = xGradient.Pixel(x,y,0) * yGradient.Pixel(x,y,0);
			H_c.Pixel(x,y,0) = yGradient.Pixel(x,y,0) * yGradient.Pixel(x,y,0);
        }
    }

    Convolve(H_a, H_a, gaussianMask);
    Convolve(H_b, H_b, gaussianMask);
    Convolve(H_c, H_c, gaussianMask);


    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {   
			// TODO:  Compute the harris score for 'srcImage' at this pixel and store in 'harrisImage'.  See the project
            //   page for pointers on how to do this
            //double determinantH = H_a.Pixel(x,y,0) * H_c.Pixel(x,y,0) - H_b.Pixel(x,y,0)*H_b.Pixel(x,y,0);
            //double traceH = H_a.Pixel(x,y,0) + H_c.Pixel(x,y,0);
            harrisImage.Pixel(x,y,0) = (H_a.Pixel(x,y,0) * H_c.Pixel(x,y,0) - H_b.Pixel(x,y,0)*H_b.Pixel(x,y,0))/ ( H_a.Pixel(x,y,0) + H_c.Pixel(x,y,0));
            
        }
    }
}



// TO DO---------------------------------------------------------------------
// Loop through the harrisImage to threshold and compute the local maxima in a neighborhood
// srcImage:  image with Harris values
// destImage: Assign 1 to a pixel if it is above a threshold and is the local maximum in 3x3 window, 0 otherwise.
//    You'll need to find a good threshold to use.
void computeLocalMaxima(CFloatImage &srcImage,CByteImage &destImage)
{
	float threshold = 0.05;

	int w = srcImage.Shape().width;
	int h = srcImage.Shape().height;

	for (int y = 0; y < h; y++) 
	{
        for (int x = 0; x < w; x++) 
		{
			float center = srcImage.Pixel(x,y,0);
			if(center >= threshold)
				destImage.Pixel(x,y,0) = 1;
			else
				destImage.Pixel(x,y,0) = 0;
			
			for(int i = 0; i < 5; i++)
			{
				for(int j = 0; j < 5; j++)
				{
					int xx = x + i - 2;
					int yy = y + j - 2;
					float value;
					if(yy == y && xx == x)
					{
						continue;
					}
					//edge cases
					if(yy < 0 || xx < 0 || yy >= h || xx >= w)
					{
						value = center;
					}
					else
					{
		 				value = srcImage.Pixel(xx, yy, 0);
					}
					if(center < value)
						destImage.Pixel(x,y,0)=0;
				}
			}
		}
	}
}
//Compute the descriptor
//The very simple 5*5 window decriptor
void computeSimpleDescriptors(CFloatImage &grayImage, FeatureSet &features) {
	vector<Feature>::iterator i = features.begin();
	while (i != features.end()) {
		Feature &f = *i;
		int w = grayImage.Shape().width;
		int h = grayImage.Shape().height;

		f.data.resize(5*5);

		for (int j = -2; j < 3; j++) {
			for (int k = -2; k < 3; k++) {
				if (f.x+k >= 0 && f.x+k < w && f.y+j >= 0 && f.y+j < h) {
					f.data.push_back(grayImage.Pixel(f.x+k, f.y+j, 1));
				} else {
					f.data.push_back(0);
				}
			}
		}

		i++;

	}
}

void computeMOPsDescriptors(CFloatImage &grayImage, FeatureSet &features) {
	
	

	int w = grayImage.Shape().width;
    int h = grayImage.Shape().height;
	CFloatImage sampleImage(8,8,1);

	// //Blur with Gaussian mask
	CFloatImage blurImage(w, h,1);
	CFloatImage gaussianMask(5,5,1);
 
    for (int i = 0; i < 5; i++)
    {
        for(int j = 0; j < 5; j++)
        {
            gaussianMask.Pixel(i,j,0) = gaussian5x5[i * 5 + j];
        }
    }
    printf("Things going well\n");
    Convolve(grayImage, blurImage, gaussianMask);

    

    // Compute the gradient (for orientation) 
	
	CFloatImage xGradient(grayImage.Shape());
	CFloatImage yGradient(grayImage.Shape());

	Convolve(grayImage, xGradient, ConvolveKernel_SobelX);
	Convolve(grayImage, yGradient, ConvolveKernel_SobelY);
    
    

    vector<Feature>::iterator i = features.begin();
 //    //TEST FOR DEBUG
	// while (i != features.end()) {
	// 	Feature &f = *i;
	// 	int w = grayImage.Shape().width;
	// 	int h = grayImage.Shape().height;

	// 	f.data.resize(5*5);

	// 	for (int j = -2; j < 3; j++) {
	// 		for (int k = -2; k < 3; k++) {
	// 			if (f.x+k >= 0 && f.x+k < w && f.y+j >= 0 && f.y+j < h) {
	// 				f.data.push_back(grayImage.Pixel(f.x+k, f.y+j, 1));
	// 			} else {
	// 				f.data.push_back(0);
	// 			}
	// 		}
	// 	}

	// 	i++;

	// }



	while (i != features.end()) {
		Feature &f = *i;

		//Compute the orientation
		f.angleRadians=atan(yGradient.Pixel(f.x,f.y,0)/xGradient.Pixel(f.x,f.y,0));
			

		CTransform3x3 xform;
		
		CTransform3x3 origin;
		origin = origin.Translation(float (-8/2), float(-8/2));

		CTransform3x3 rotate;
        rotate = rotate.Rotation(f.angleRadians * 180.0 / PI);

        CTransform3x3 trans;
        trans = trans.Translation( f.x, f.y );

        //using spacing of s=5 pixels between samples
        CTransform3x3 scale;
        scale[0][0] = 5;
        scale[1][1] = 5;

        xform = trans * rotate * scale * origin;

        //rotate to horizontal
        WarpGlobal(blurImage, sampleImage, xform ,eWarpInterpLinear);

        f.data.resize(8*8);


        //After sampling, the descriptor vector is nomalised
        double mean = 0.0;
        for (int row = 0; row < 8; row++ ){
            for (int col = 0; col < 8; col++){
                mean += sampleImage.Pixel(col, row, 0);
            }
        }
        mean = mean / (8 * 8);


        double variance = 0.0;
        for (int row = 0; row < 8; row++ ){
            for (int col = 0; col < 8; col++){
                variance += pow((sampleImage.Pixel(col, row, 0) - mean), 2);
            }
        }
        double stndiv = sqrt( variance / ( (8 * 8) - 1));


        //Fill in the data
        for (int row = 0; row < 8; row++ ){
            for (int col = 0; col < 8; col++){
                f.data[row * 8 + col] = (sampleImage.Pixel(col, row, 0) - mean)/stndiv;
            }
        }


		i++;
	}
	

}



// Perform simple feature matching.  This just uses the SSD
// distance between two feature vectors, and matches a feature in the
// first image with the closest feature in the second image.  It can
// match multiple features in the first image to the same feature in
// the second image.
void ssdMatchFeatures(const FeatureSet &f1, const FeatureSet &f2, vector<FeatureMatch> &matches, double &totalScore) {
	int m = f1.size();
	int n = f2.size();

	matches.resize(m);
	totalScore = 0;

	double d;
	double dBest;
	int idBest;

	for (int i=0; i<m; i++) {
		dBest = 1e100;
		idBest = 0;

		for (int j=0; j<n; j++) {
			d = distanceSSD(f1[i].data, f2[j].data);

			if (d < dBest) {
				dBest = d;
				idBest = f2[j].id;
			}
		}

        matches[i].id1 = f1[i].id;
		matches[i].id2 = idBest;
		matches[i].score = dBest;
		totalScore += matches[i].score;
	}
}

// TODO: Write this function to perform ratio feature matching.  
// This just uses the ratio of the SSD distance of the two best matches as the score
// and matches a feature in the first image with the closest feature in the second image.
// It can match multiple features in the first image to the same feature in
// the second image.  (See class notes for more information, and the sshMatchFeatures function above as a reference)
void ratioMatchFeatures(const FeatureSet &f1, const FeatureSet &f2, vector<FeatureMatch> &matches, double &totalScore) 
{
	int m = f1.size();
	int n = f2.size();

	matches.resize(m);
	totalScore = 0;

	double d;
	double dBest;
	double secondBest;
	int idBest;
	int idSecondbest;

     
	for (int i = 0; i < m; i++) {
		dBest = 1e100;
		idBest = 0;

		for (int j = 0; j < n; j++) {
			d = distanceSSD(f1[i].data, f2[j].data);

			if (d < dBest) {
				secondBest = dBest;
				idSecondbest = idBest;
				dBest = d;
				idBest = f2[j].id;
			}
		}

        matches[i].id1 = f1[i].id;
		matches[i].id2 = idBest;
		
		matches[i].score = dBest / secondBest;

		totalScore += matches[i].score;
	}
    
}


// Convert Fl_Image to CFloatImage.
bool convertImage(const Fl_Image *image, CFloatImage &convertedImage) {
	if (image == NULL) {
		return false;
	}

	// Let's not handle indexed color images.
	if (image->count() != 1) {
		return false;
	}

	int w = image->w();
	int h = image->h();
	int d = image->d();

	// Get the image data.
	const char *const *data = image->data();

	int index = 0;

	for (int y=0; y<h; y++) {
		for (int x=0; x<w; x++) {
			if (d < 3) {
				// If there are fewer than 3 channels, just use the
				// first one for all colors.
				convertedImage.Pixel(x,y,0) = ((uchar) data[0][index]) / 255.0f;
				convertedImage.Pixel(x,y,1) = ((uchar) data[0][index]) / 255.0f;
				convertedImage.Pixel(x,y,2) = ((uchar) data[0][index]) / 255.0f;
			}
			else {
				// Otherwise, use the first 3.
				convertedImage.Pixel(x,y,0) = ((uchar) data[0][index]) / 255.0f;
				convertedImage.Pixel(x,y,1) = ((uchar) data[0][index+1]) / 255.0f;
				convertedImage.Pixel(x,y,2) = ((uchar) data[0][index+2]) / 255.0f;
			}

			index += d;
		}
	}
	
	return true;
}

// Convert CFloatImage to CByteImage.
void convertToByteImage(CFloatImage &floatImage, CByteImage &byteImage) {
	CShape sh = floatImage.Shape();

    assert(floatImage.Shape().nBands == byteImage.Shape().nBands);
	for (int y=0; y<sh.height; y++) {
		for (int x=0; x<sh.width; x++) {
			for (int c=0; c<sh.nBands; c++) {
				float value = floor(255*floatImage.Pixel(x,y,c) + 0.5f);

				if (value < byteImage.MinVal()) {
					value = byteImage.MinVal();
				}
				else if (value > byteImage.MaxVal()) {
					value = byteImage.MaxVal();
				}

				// We have to flip the image and reverse the color
				// channels to get it to come out right.  How silly!
				byteImage.Pixel(x,sh.height-y-1,sh.nBands-c-1) = (uchar) value;
			}
		}
	}
}

// Compute SSD distance between two vectors.
double distanceSSD(const vector<double> &v1, const vector<double> &v2) {
	int m = v1.size();
	int n = v2.size();

	if (m != n) {
		// Here's a big number.
		return 1e100;
	}

	double dist = 0;

	for (int i=0; i<m; i++) {
		dist += pow(v1[i]-v2[i], 2);
	}

	
	return sqrt(dist);
}

// Transform point by homography.
void applyHomography(double x, double y, double &xNew, double &yNew, double h[9]) {
	double d = h[6]*x + h[7]*y + h[8];

	xNew = (h[0]*x + h[1]*y + h[2]) / d;
	yNew = (h[3]*x + h[4]*y + h[5]) / d;
}

// Compute AUC given a ROC curve
double computeAUC(vector<ROCPoint> &results)
{
	double auc=0;
	double xdiff,ydiff;
	for (int i = 1; i < (int) results.size(); i++)
    {
        //fprintf(stream,"%lf\t%lf\t%lf\n",thresholdList[i],results[i].falseRate,results[i].trueRate);
		xdiff=(results[i].falseRate-results[i-1].falseRate);
		ydiff=(results[i].trueRate-results[i-1].trueRate);
		auc=auc+xdiff*results[i-1].trueRate+xdiff*ydiff/2;
    	    
    }
	return auc;
}