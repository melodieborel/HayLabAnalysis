import deeplabcut

path_config_file="/mnt/data/AurelieB_dlc/CheeseboardMiniscope-AurelieB-2025-02-21/config.yaml"

#deeplabcut.create_training_dataset(path_config_file,Shuffles=[1])
#deeplabcut.train_network(path_config_file, shuffle=1, saveiters=1000, displayiters=10)

deeplabcut.evaluate_network(path_config_file, plotting=True)