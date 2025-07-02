import os
import pandas as pd

# Get the path to the directory containing the files
data_dir = 'xxxxxxxx/crop_part1'

# Get a list of all the files in the directory
files = os.listdir(data_dir)

# Create empty lists to store the age, gender, and race data
age_list = []
gender_list = []
race_list = []

# Loop through each file in the directory
for file in files:
  # Split the file name into parts using the '_' character
  parts = file.split('_')

  # Check if the file name has at least three '_'s
  if len(parts) >= 3:

    try:
        age = int(parts[0])
        gender = int(parts[1])
        race = int(parts[2])

        # Add the extracted information to the lists
        age_list.append(age)
        gender_list.append(gender)
        race_list.append(race)
    except:
        print(file)
  else:
    print(file)

# Create a Pandas DataFrame with the age, gender, and race data
df = pd.DataFrame({
  'age': age_list,
  'gender': gender_list,
  'race': race_list
})

bins = [0, 20, 40, 60, 80, 116]

# Create labels for age groups
labels = ['0-20', '21-40', '41-60', '61-80', '81-116']

# Convert age to age groups
df['age_group'] = pd.cut(df['age'], bins, labels=labels)

# print(df.groupby('age_group').size().to_frame(name='count'))
# print(df.groupby('gender').size().to_frame(name='count'))
# print(df.groupby('race').size().to_frame(name='count'))

pivot_table1 = df.pivot_table(index='age_group', columns=['gender'], aggfunc=len).transpose()
# Print the pivot table
print(pivot_table1)
pivot_table2 = df.pivot_table(index='age_group', columns=['race'], aggfunc=len).transpose()
print(pivot_table2)
pivot_table1.to_excel('xxxxx.xlsx')
pivot_table2.to_excel('xxxxxxx2.xlsx')
