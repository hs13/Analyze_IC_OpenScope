% correlation between SVM score and neural activity

trialsoi = SVMpsth.trialorder==106;
figure; plot(psthbintli, squeeze(mean(SVMpsth.score(:,trialsoi,1,:),2)))

figure; hold all
for typi = 1:4
plot(psthbintli, squeeze(mean(SVMpsth.score(:,trialsoi,typi,:),[2,4])))
end
