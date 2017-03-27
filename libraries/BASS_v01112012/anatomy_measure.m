function data_out = anatomy_measure(data,region_vert,op,reshape_mat)


[Nr Nv] = size(region_vert);

data_out = zeros(size(data));
data_out(Nr+1:end,:) = [];

if nargin<4
   reshape_mat = size(data_out); 
end

switch op
    case 'mean'
        
        for ir=1:Nr
            data_out(ir,:) = mean( data(region_vert(ir,:),:) ,1);
        end
        
    case 'sum'
        
        disp('Not implmented yet - sorry')
        data_out = [];
        
    case 'zscore'
        
        disp('Not implmented yet - sorry')
        data_out = [];
        
    otherwise
        error('Not an option.')
end

data_out = reshape(data_out,reshape_mat);