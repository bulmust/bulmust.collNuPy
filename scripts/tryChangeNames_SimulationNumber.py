# import OS module
import os as osCommand
import shutil

def listOfFoldersInDirectory(path):
    return [f for f in osCommand.listdir(path) \
        if osCommand.path.isdir(osCommand.path.join(path, f))]
if __name__ == '__main__':
    # Current directory
    currentDir = osCommand.getcwd()
    for it_sim in listOfFoldersInDirectory(currentDir):
        if 'E6_11' in it_sim: # Second 30 Simulation
            resultsFolder= it_sim+'/results'
            resultsFolder_inside = listOfFoldersInDirectory(resultsFolder)
            # Sort the list
            resultsFolder_inside.sort(key=lambda x: int(x.split('simulation')[-1]))
            
            # Remove simulation1, simulation2 and simulation3 folders
            for it in resultsFolder_inside[0:3]:
                shutil.rmtree(resultsFolder+'/'+it)
            # Remove simulation1, simulation2 and simulation3 from the list
            resultsFolder_inside = resultsFolder_inside[3:]

            for it_sim_inside in resultsFolder_inside:
                # Rename the folder
                osCommand.rename(resultsFolder+'/'+it_sim_inside,resultsFolder+'/simulation'+str(int(it_sim_inside.split('simulation')[-1])+87))
        elif 'E11_17' in it_sim: # Third 30 Simulation
            resultsFolder= it_sim+'/results'
            resultsFolder_inside = listOfFoldersInDirectory(resultsFolder)
            # Sort the list
            resultsFolder_inside.sort(key=lambda x: int(x.split('simulation')[-1]))
            
            # Remove simulation1, simulation2 and simulation3 folders
            for it in resultsFolder_inside[0:3]:
                shutil.rmtree(resultsFolder+'/'+it)
            # Remove simulation1, simulation2 and simulation3 from the list
            resultsFolder_inside = resultsFolder_inside[3:]

            for it_sim_inside in resultsFolder_inside:
                # Rename the folder
                osCommand.rename(resultsFolder+'/'+it_sim_inside,resultsFolder+'/simulation'+str(int(it_sim_inside.split('simulation')[-1])+174))
        elif 'E17_23' in it_sim: # Fourth 30 Simulation
            resultsFolder= it_sim+'/results'
            resultsFolder_inside = listOfFoldersInDirectory(resultsFolder)
            # Sort the list
            resultsFolder_inside.sort(key=lambda x: int(x.split('simulation')[-1]))
            # Remove simulation1, simulation2 and simulation3 folders
            for it in resultsFolder_inside[0:3]:
                shutil.rmtree(resultsFolder+'/'+it)
            # Remove simulation1, simulation2 and simulation3 from the list
            resultsFolder_inside = resultsFolder_inside[3:]

            for it_sim_inside in resultsFolder_inside:
                # Rename the folder
                osCommand.rename(resultsFolder+'/'+it_sim_inside,resultsFolder+'/simulation'+str(int(it_sim_inside.split('simulation')[-1])+261))
        elif 'E23_29' in it_sim: # Fifth 30 Simulation
            resultsFolder= it_sim+'/results'
            resultsFolder_inside = listOfFoldersInDirectory(resultsFolder)
            # Sort the list
            resultsFolder_inside.sort(key=lambda x: int(x.split('simulation')[-1]))
            
            # Remove simulation1, simulation2 and simulation3 folders
            for it in resultsFolder_inside[0:3]:
                shutil.rmtree(resultsFolder+'/'+it)
            # Remove simulation1, simulation2 and simulation3 from the list
            resultsFolder_inside = resultsFolder_inside[3:]

            for it_sim_inside in resultsFolder_inside:
                # Rename the folder
                osCommand.rename(resultsFolder+'/'+it_sim_inside,resultsFolder+'/simulation'+str(int(it_sim_inside.split('simulation')[-1])+348))
        elif 'E29_35' in it_sim: # Sixth 30 Simulation
            resultsFolder= it_sim+'/results'
            resultsFolder_inside = listOfFoldersInDirectory(resultsFolder)
            # Sort the list
            resultsFolder_inside.sort(key=lambda x: int(x.split('simulation')[-1]))
            
            # Remove simulation1, simulation2 and simulation3 folders
            for it in resultsFolder_inside[0:3]:
                shutil.rmtree(resultsFolder+'/'+it)
            # Remove simulation1, simulation2 and simulation3 from the list
            resultsFolder_inside = resultsFolder_inside[3:]

            for it_sim_inside in resultsFolder_inside:
                # Rename the folder
                osCommand.rename(resultsFolder+'/'+it_sim_inside,resultsFolder+'/simulation'+str(int(it_sim_inside.split('simulation')[-1])+435))
        elif 'E35_41' in it_sim: # Seventh 30 Simulation
            resultsFolder= it_sim+'/results'
            resultsFolder_inside = listOfFoldersInDirectory(resultsFolder)
            # Sort the list
            resultsFolder_inside.sort(key=lambda x: int(x.split('simulation')[-1]))
            
            # Remove simulation1, simulation2 and simulation3 folders
            for it in resultsFolder_inside[0:3]:
                shutil.rmtree(resultsFolder+'/'+it)
            # Remove simulation1, simulation2 and simulation3 from the list
            resultsFolder_inside = resultsFolder_inside[3:]

            for it_sim_inside in resultsFolder_inside:
                # Rename the folder
                osCommand.rename(resultsFolder+'/'+it_sim_inside,resultsFolder+'/simulation'+str(int(it_sim_inside.split('simulation')[-1])+522))
        