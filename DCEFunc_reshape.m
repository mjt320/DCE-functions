function dataOut=DCEFunc_reshape(dataIn,dims)
%convert a 4D array to a 2D array, such that each column represents
%a time series...or the reverse of the above. Useful for converting a 4D
%image into a many time series for processing.
%INPUT:
%dataIn: 4D image is converted into a matrix (dataOut), in which each row
%corresponds to a time point and each column corresponds to a 3D voxel.
%
%Alternatively a time x voxel 2D array (dataIn) is converted back to a 4D
%image (dataOut) using the sizes specified for the first 3 dimensions (dims). To reverse
%the 4D -> 2D operation, dims should be a vector containing the first 3
%dimensions of the 4D image.
%
%For example, DCEFunc_reshape(DCEFunc_reshape(a4D),[size(a4D,1) size(a4D,2)
%size(a4D,3)]) should return a4D

switch ndims(dataIn)
    case 2
        NFrames=size(dataIn,1);
        dataOut=nan([dims NFrames]);
        for iFrame=1:NFrames
            dataOut(:,:,:,iFrame)=reshape(dataIn(iFrame,:),dims);
        end
    case 4
        NFrames=size(dataIn,4);
        dataOut=nan(NFrames,size(dataIn,1)*size(dataIn,2)*size(dataIn,3));
        for iFrame=1:NFrames
            temp=dataIn(:,:,:,iFrame);
            dataOut(iFrame,:)=temp(:);
        end
    otherwise
        error('Input must be 2 or 4 dimensional.');
end
