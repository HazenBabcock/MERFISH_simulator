#!/usr/bin/env python
"""
Creates the movies for each FOV and imaging cycle.
"""
import os


def createMovie(config, simParams, dataPath, fov, iRound):
    """
    Creates the movie for a specific fov and imaging round.
    """
    do = simParams.get_data_organization()
    pos = simParams.get_fov_xy_um(fov)

    # Create images for movie.
    #
    # This first determines the frame type for each frame in
    # the movie based on the data organization file. Then it
    # creates the image using simulation tasks that end in
    # 'image'. At least one of these much match the 'channelName'
    # field for this frame in the data organization file.
    #
    images = []
    for fi in range(do.get_number_frames(iRound)):
        desc = do.get_frame_description(iRound, fi)

        image = None
        for elt in config:
            if elt.endswith("image"):
                tmp = config[elt].make_image(config,
                                             simParams,
                                             fov,
                                             iRound,
                                             desc)
                if tmp is not None:
                    if image is None:
                        image = tmp
                    else:
                        image += tmp

        assert (image is not None), "No image created for " + desc[0] + "!"
        
        images.append(image)

    # Process with camera.
    #
    stack = []
    for elt in images:
        stack.append(config["camera"].camera_image(elt))

    # Save.
    #
    # FIXME: This should use the imageRegExp to figure out the name.
    #
    if (iRound >= 0):
        fileName = do.get_image_type(iRound) + "_{0:03d}_{1:02d}".format(fov, iRound)
    else:
        fileName = do.get_image_type(iRound) + "_{0:03d}".format(fov)
        
    filePath = os.path.join(dataPath, fileName)    
    config["movie_writer"].save_stack(filePath, stack, pos)
