
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: eigfaces.cpp                                                                         //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////



#include "stdafx.h"
#include <algorithm>


EigFaces::EigFaces()
:
Faces()
{
	//empty
}

EigFaces::EigFaces(int count, int width, int height)
:
Faces(count, width, height)
{
	//empty
}

void EigFaces::projectFace(const Face& face, Vector& coefficients) const
{
	if (face.getWidth()!=width || face.getHeight()!=height) {
		throw Error("Project: Face to project has different dimensions");
	}

	coefficients.resize(getSize());
	// ----------- TODO #2: compute the coefficients for the face and store in coefficients.
	Face difference = Face(width, height);
	face.sub(average_face, difference);
	
	//std::cout << average_face.getSize() << std::endl;
	//std::cout<< difference.getSize() << std::endl;

	for (int i = 0; i < getSize(); i++) {
		coefficients[i] = difference.dot((*this)[i]);
	}


}

void EigFaces::constructFace(const Vector& coefficients, Face& result) const
{	
	// ----------- TODO #3: construct a face given the coefficients
	//x = miu + u1*w1+u2*w2+...
	for (int j = 0; j < result.getSize(); ++j)
	{
		result[j] = average_face[j];//miu
		for (int i = 0; i < getSize(); ++i)
		{
			result[j] += (*this)[i][j] * coefficients[i];
		}
	}

}

bool EigFaces::isFace(const Face& face, double max_reconstructed_mse, double& mse) const
{
	// ----------- TODO #4: Determine if an image is a face and return true if it is. Return the actual
	// MSE you calculated for the determination in mse
	// Be sure to test this method out with some face images and some non face images
	// to verify it is working correctly.
	Vector coefficients;
	Face constructedFace = Face(width, height);
	projectFace(face, coefficients);
	constructFace(coefficients, constructedFace);
	mse = face.mse(constructedFace);
	
	if (mse < max_reconstructed_mse) {
		return true;
	} else {
		return false;
	}
}

bool EigFaces::verifyFace(const Face& face, const Vector& user_coefficients, double max_coefficients_mse, double& mse) const
{
	// ----------- TODO #5 : Determine if face is the same user give the user's coefficients.
	// return the MSE you calculated for the determination in mse.
	Vector calCoeff;
	projectFace(face, calCoeff);
	mse = calCoeff.mse(user_coefficients);
	if (mse < max_coefficients_mse) {
		return true;
	} else {
		return false;
	}
	
}

void EigFaces::recognizeFace(const Face& face, Users& users) const
{
	// ----------- TODO #6: Sort the users by closeness of match to the face
	Vector coefficients;
	projectFace(face,coefficients);
	for (int i = 0; i < users.getSize(); i++) {
		double mse = users[i].mse(coefficients);
		
		users[i].setMse(mse);
	}
	users.sort();
}

bool compError(const FacePosition & a, const FacePosition & b) { return a.error < b.error; }

void EigFaces::findFace(const Image& img, double min_scale, double max_scale, double step, int n, bool crop, Image& result) const
{
	// ----------- TODO #7: Find the faces in Image. Search image scales from min_scale to max_scale inclusive,
	// stepping by step in between. Find the best n faces that do not overlap each other. If crop is true,
	// n is one and you should return the cropped original img in result. The result must be identical
	// to the original besides being cropped. It cannot be scaled and it must be full color. If crop is
	// false, draw green boxes (use r=100, g=255, b=100) around the n faces found. The result must be
	// identical to the original image except for the addition of the boxes.
	
	
	
	Face  test(width, height);

	std::list<FacePosition> bestPos;
	std::list<FacePosition>::iterator pos;


	for (double scale = min_scale; scale <= max_scale; scale += step) {
		Image scaledImg(scale * img.getWidth(), scale * img.getHeight(), img.getColors()); 
		img.resample(scaledImg);

		for (int x = 0; x < scaledImg.getWidth() - width; x++) {
			for (int y = 0; y < scaledImg.getHeight() - height; y++) {
				test.subimage(x, x + test.getWidth()-1, y, y + test.getHeight()-1, scaledImg, false);

				double mse;
				bool isface = isFace(test, 1000.0, mse);

				//If the parameter crop is true, then you should return the image cropped to the best face in result. 
				//Note that this result must be the same scale and resolution as the original image. 
				//If crop is false, then draw bright green boxes around each found face in the original image and return this in result.

				//need to find face first
				//Whenever a position is found that is better than the last position in the list, 
				//insert it into the proper place in the list such that the list stays sorted, 
				//and remove the last element in the list.

				if (isface) {
					FacePosition currentPos = FacePosition();
					currentPos.x = x;
					currentPos.y = y;
					currentPos.scale = scale;
					currentPos.error = mse;

					pos = bestPos.begin();

					
					if (bestPos.empty()) {
						bestPos.push_back(currentPos);
					} else {
						//if there is already a face position in the list that overlaps the one you just found. 
						//In that case you must choose the one with the better face and throw the other away.
						

						bool overlap = false;
						//iterate and check whether overlap
						for(pos = bestPos.begin(); pos != bestPos.end(); pos++) {
							if (abs((*pos).x/(*pos).scale - x/scale) < width/(*pos).scale
							 || abs((*pos).y/(*pos).scale - y/scale) < height/(*pos).scale
							 || abs((*pos).x/(*pos).scale - x/scale) < width/scale
							 || abs((*pos).y/(*pos).scale - y/scale) < height/scale) {
								//choose the one with the better face
								overlap = true;
								if(currentPos.error < (*pos).error){
									(*pos).x = x;
									(*pos).y = y;
									(*pos).scale = scale;
									(*pos).error = mse;
									//*pos = currentPos;
									
									bestPos.sort(compError);
								}
								break;
							}
						}
						
						if (!overlap) {
							//check whether full
							//if (bestPos.size() < n) {
								bestPos.push_back(currentPos);
								bestPos.sort(compError);
							//}
							// else {
							// 	//if full, compare with the last, and replace if currentPos is better
							// 	FacePosition worst = *bestPos.end();
							// 	if (currentPos.error < worst.error) {
							// 		bestPos.pop_back();
							// 		bestPos.push_back(currentPos);
							// 		bestPos.sort(compError);
							// 	}

							// }
						}

					}
				}

				
			}
		}
	}

	while(bestPos.size()>n){
		bestPos.pop_back();
	}
		
	if (crop) {
		FacePosition best = bestPos.front();
			
		std::cout << "scale: " << best.scale << std::endl;
			
		img.crop(best.x / best.scale, best.y/best.scale,(best.x + width) / best.scale, (best.y + height)/best.scale, result);

	} else {
		result.resize(img.getWidth(),img.getHeight(),img.getColors());
		img.resample(result);

		int min_x,min_y,max_x,max_y;

		for(pos = bestPos.begin(); pos != bestPos.end(); pos++) {
			min_x = (int)((*pos).x / (*pos).scale);
			min_y = (int)((*pos).x / (*pos).scale);
			max_x = (int)(((*pos).x + width) / (*pos).scale);
			max_y = (int)(((*pos).y + width*10/7) / (*pos).scale);

			result.line(min_x, min_y, max_x, min_y, 100, 255, 100);
			result.line(min_x, max_y, max_x, max_y, 100, 255, 100);
			result.line(min_x, min_y, min_x, max_y, 100, 255, 100);
			result.line(max_x, min_y, max_x, max_y, 100, 255, 100);
		}

	}

	

}

void EigFaces::morphFaces(const Face& face1, const Face& face2, double distance, Face& result) const
{
	// TODO (extra credit): MORPH along *distance* fraction of the vector from face1 to face2 by
	// interpolating between the coefficients for the two faces and reconstructing the result.
	// For example, distance 0.0 will approximate the first, while distance 1.0 will approximate the second.
	// Negative distances are ok two.

}

const Face& EigFaces::getAverage() const
{
	return average_face;
}

void EigFaces::setAverage(const Face& average)
{
	average_face=average;
}



