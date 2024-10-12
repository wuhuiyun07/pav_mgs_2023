import os

# Specify the file path

parent_directory = 'results/megahit/'
SAMPLES = [f.name for f in os.scandir(parent_directory) if f.is_dir()]
subfolder_count = sum(1 for f in os.scandir(parent_directory) if f.is_dir()) ##78 subfolders, samples

print(SAMPLES)
print(subfolder_count)
# ['Wu_18_1_S6', 'Wu_29_2_S22', 'Wu_16_1_S1', 'Wu_17_3_S3', 'Wu_17_4_S4', 'Wu_16_2_S2', 'Wu_8_4_S10', 'Wu_24_5_S20', 'Wu_14_3_S23', 'Wu_23_1_S11', 'Wu_21_1_S16', 'Wu_19_3_S13', 'Wu_23_5_S15', 'Wu_29_5_S25', 'Wu_29_3_S23', 'Wu_21_2_S17', 'Wu_19_5_S15', 'Wu_21_5_S20', 'Wu_29_1_S21', 'Wu_9_3_S14', 'Wu_10_3_S19', 'Wu_4_1_S2', 'Wu_26_3_S23', 'Wu_4_3_S4', 'Wu_19_1_S11', 'Wu_24_2_S17', 'Wu_14_1_S21', 'Wu_18_3_S8', 'Wu_9_4_S15', 'Wu_23_2_S12', 'Wu_17_5_S5', 'Wu_24_1_S16', 'Wu_9_2_S13', 'Wu_10_2_S18', 'Wu_8_5_S11', 'Wu_26_4_S24', 'Wu_22_4_S9', 'Wu_8_2_S8', 'Wu_NEG_S1', 'Wu_14_4_S24', 'Wu_24_3_S18', 'Wu_21_3_S18', 'Wu_4_4_S5', 'Wu_14_5_S25', 'Wu_16_4_S4', 'Wu_8_1_S7', 'Wu_22_3_S8', 'Wu_22_1_S6', 'Wu_24_4_S19', 'Wu_16_5_S5', 'Wu_19_2_S12', 'Wu_26_2_S22', 'Wu_17_2_S2', 'Wu_8_3_S9', 'Wu_9_5_S16', 'Wu_18_5_S10', 'Wu_22_5_S10', 'Wu_23_3_S13', 'Wu_9_1_S12', 'Wu_14_2_S22', 'Wu_4_5_S6', 'Wu_18_2_S7', 'Wu_29_4_S24', 'Wu_26_5_S25', 'Wu_4_2_S3', 'Wu_16_3_S3', 'Wu_21_4_S19', 'Wu_23_4_S14', 'Wu_19_4_S14', 'Wu_Blank032824_S27', 'Wu_17_1_S1', 'Wu_10_1_S17', 'Wu_18_4_S9', 'Wu_Positive_S26', 'Wu_10_4_S20', 'Wu_26_1_S21', 'Wu_22_2_S7', 'Wu_POS_S26']