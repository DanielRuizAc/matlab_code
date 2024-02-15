load("DataSets\OnlyXsens\transTrainData.mat");

%% Create complete Dataset

without_exo_keys = without_exo_data.keys;
Comp_dataset = [];
for key_id=1:length(without_exo_keys)
    key = without_exo_keys(key_id);
    Comp_dataset = [Comp_dataset; without_exo_data(key).filtered_db];
end

with_exo_keys = with_exo_data.keys;
for key_id=1:length(with_exo_keys)
    key = with_exo_keys(key_id);
    Comp_dataset = [Comp_dataset; with_exo_data(key).filtered_db];
end
%% Normalize and apply pca

[dataset_norm, C, S] = normalize(Comp_dataset, "zscore");
[coeff_pca, score,latent] = pca(dataset_norm);

dataset_pca = dataset_norm;
for i=1:size(dataset_norm,1)
    dataset_pca(i,:) = coeff_pca'*(dataset_norm(i,:)');
end

%% Applying KMeans

[ids_eucl, C_eucl] = kmeans(dataset_pca(:,1),2,"Distance","sqeuclidean");
[ids_city, C_city] = kmeans(dataset_pca(:,1),2,"Distance","cityblock");
[ids_cos, C_cos] = kmeans(dataset_pca(:,1),2,"Distance","cosine");

[ids_medoids_sqeucl, C_medoids_sqeucl] = kmedoids(dataset_pca(:,1),2,"Distance","sqeuclidean");
% ids_spectral = spectralcluster(dataset_pca(:,1),2);
C_eucl = sort(C_eucl);
% save(".\DataSets\OnlyXsens\kmeans_data", "C", "S", "coeff_pca", "C_eucl", "C_city", "C_cos", "C_medoids_sqeucl")

%% Plot results with training data
figure
plot(ids_eucl)
% hold on
% plot(ids_city)
% hold on
% plot(ids_cos)
% hold on
% plot(ids_medoids_eucl)
% legend
% hold on
% plot(ids_eucl)