We need a way to exchange information between the matlab process running the full T&C model, potentially online, and the python process, training and inferring using the stomatal resistance model.

# Using sockets
Matlab apparently has an interface for TCP communication [here](https://www.mathworks.com/help/instrument/tcp-ip-interface.html). That could be a reasonnably sound solution, but would require coming up with a protocol, and might require a significant synchronisation effort.

# Using files
Both process could write partial results on specific files. That would be a relatively inefficient solution, and would require some sort of polling mechanism to notice updates? Apparently, file watchers exist in python.

# Using matlab sessions / python matlab-engine
There seems to be an existing solution to combine both languages, presented [here](https://www.mathworks.com/help/matlab/matlab_external/connect-python-to-running-matlab-session.html). This might be the cleanest solution, but would require a better command of matlab on my end. Also, I don't know if running a very large process like T&C would fit this workflow well. [this](https://blogs.mathworks.com/student-lounge/2020/09/14/using-matlab-and-python-together/) also seems like a good resource. Presents both way (calling MATLAB from python and vice-versa). It would probably require converting the python code to a library format. But then, it could become hard to have something like a FCN model loaded in memory, if we have to call a purely functional function every time. [here's](https://www.mathworks.com/videos/mathworks-energy-conference-2022-how-python-and-matlab-can-work-together-1679461122169.html) another resource, video.


The solution we choose will depend on both the latency and throughput requirements. If those are low, a very simple solution will be enough.

