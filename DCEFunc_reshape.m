function dataOut=DCEFunc_reshape(dataIn,dims)
%convert a 4D array to a 2D array, such that each column represents
%a time series
%...or the reverse of the above
%for the reverse operation, the row vector argument dim must be used defining the size of the first
%three dimensions of the 4D output

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
