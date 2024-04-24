
vect = [randi([0,1],5,1)';
    randi([0,1],5,1)';
    randi([0,1],5,1)';]; % random vector of bits



H =[1 1 0 1 1 0 0 1 0 0;
    0 1 1 0 1 1 1 0 0 0;
    0 0 0 1 0 0 0 1 1 1;
    1 1 0 0 0 1 1 0 1 0;
    0 0 1 0 0 1 0 1 0 1];

% H = rref(H)
% mod(H,2)

P =[0     1     1     1     0;
    1     0     1     0     0;
    1     0     1     0     1;
    0     0     1     1     1;
    1     1     0     0     1].'; 

I_5 = eye(5);

Gen = [P I_5];

vect_enc = vect * Gen;

vect_enc = reshape(vect_enc',[1 size(vect_enc,1) * size(vect_enc,2)]);
for m = 1:length(vect_enc)
    vect_enc(m) = mod(vect_enc(m),2);
    
end
vect_enc_with_error = vect_enc;
vect_enc_with_error(25) = vect_enc_with_error(25)+1; % insert error
vect_enc_with_error = mod(vect_enc_with_error,2);
decoded = decoding(vect_enc_with_error,P,H);



verif = decoded- vect_enc;


