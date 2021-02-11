import numpy as np
import plotly.graph_objs as go

facies_code = 0
temp_01 = np.linspace(prob[0].T[facies_code],prob[1].T[facies_code],50)
temp_02 = np.linspace(prob[1].T[facies_code],prob[2].T[facies_code],50)
prob_fac00 = np.append(temp_01, temp_02).reshape(100,100,100)

# Here is the code for 3d facies interpolation for facies code 1 (Shoreface):
facies_code = 1
temp_01 = np.linspace(prob[0].T[facies_code],prob[1].T[facies_code],50)
temp_02 = np.linspace(prob[1].T[facies_code],prob[2].T[facies_code],50)
prob_fac01 = np.append(temp_01, temp_02).reshape(100,100,100)


# Here is the code for 3d facies interpolation for facies code 2 (Offshore Transition):
facies_code = 2
temp_01 = np.linspace(prob[0].T[facies_code],prob[1].T[facies_code],50)
temp_02 = np.linspace(prob[1].T[facies_code],prob[2].T[facies_code],50)
prob_fac02 = np.append(temp_01, temp_02).reshape(100,100,100)

prob_stack = np.c_[prob_fac00.flatten(), prob_fac01.flatten(), prob_fac02.flatten()]

fac_transsitional = np.array([np.random.choice(a=(0,1,2),p=(i)) for i in prob_stack]).reshape(100,100,100)

grid = np.linspace(0,1,100)
grid_100 = np.linspace(0,100,99)
xx,yy = np.meshgrid(grid,grid)
xz,yz = np.meshgrid(grid,grid_100)

data_z  = go.Surface(x=xx, y=yy, z=np.full_like(xx,49), surfacecolor = fac_transsitional[49,:,:],colorscale='Viridis')
data_y  = go.Surface(x=xx, y=np.full_like(xx,0), z=yz, surfacecolor = fac_transsitional[:,0,:],colorscale='Viridis')
data_x  = go.Surface(x=np.full_like(xx,1), y=xz, z=yz, surfacecolor = fac_transsitional[:,:,99],colorscale='Viridis')
data_x1  = go.Surface(x=np.full_like(xx,.5), y=xz, z=yz, surfacecolor = fac_transsitional[:,:,49],colorscale='Viridis')

fig = go.Figure(data=[data_x, data_y, data_z, data_x1], )

fig.update_layout(
        scene = {
            "xaxis": {"nticks": 20},
            "zaxis": {"nticks": 4},
            'camera_eye': {"x": 0, "y": -0.5, "z": 0.5},
            "aspectratio": {"x": 1, "y": 1, "z": 0.3}
        },
    width=800,
    height=800,)

fig.show()
