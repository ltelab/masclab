function BER = computeBER( real, pred )
%BER Computes the BER score (balanced error rate) for a prediction given
% real labels

labels=unique(real);
n_labels=length(labels);

BER=0;
for i=1:n_labels
    nb_samples_in_class=sum(real==labels(i));
    BER=BER+1/nb_samples_in_class*sum((real==labels(i)).*(real~=pred));
end
BER=BER/n_labels;
end

