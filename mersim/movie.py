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
    images = []
    for fi in range(do.get_number_frames(iRound)):
        desc = do.get_frame_description(iRound, fi)

        if "bit" in desc[0]:
            images.append(config["bit_image"].make_image(config,
                                                         simParams,
                                                         fov,
                                                         iRound,
                                                         desc))
        else:
            task_name = desc[0] + "_image"
            images.append(config[task_name].make_image(config,
                                                       simParams,
                                                       fov,
                                                       iRound,
                                                       desc))

    # Convolve with PSF.
    stack = images


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
