figure;
ImgM=ops.meanImg;
imshow(ImgM,[])
roi_img=imshow(ImgM,[]);
imwrite(uint16(roi_img.CData),[processedDir,'mask_img_Plane.tif'])