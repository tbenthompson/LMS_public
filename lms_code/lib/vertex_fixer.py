def vertex_fixer(d):
    # Create a dictionary out of our list so that the indices of the vertices
    # don't change upon deletion
    if False:
        twod_vertices =  d['twod_vertices'].copy()
        not_done = True
        refresh = True
        del_vert_list = []
        combine_list = []
        while not_done:
            if refresh:
                refresh = False
                plt.figure(1)
                plt.clf()
                vx = [twod_vertices[key][1] for key in twod_vertices]
                vy = [twod_vertices[key][2] for key in twod_vertices]
                labels = [key for key in twod_vertices]
                plt.scatter(vx, vy, marker = 'o')
                for i, (x, y) in zip(labels, (zip(vx, vy))):
                    # Label at a radius of 20 pixels from the point in a random
                    # direction
                    rand_theta = np.random.rand() * 2 * np.pi
                    label_x_loc = cos(rand_theta) * 20
                    label_y_loc = sin(rand_theta) * 20
                    plt.annotate(i, xy = (x, y), xytext = (label_x_loc, label_y_loc),
                            textcoords = 'offset points', ha = 'right', va = 'bottom',
                            bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
                            arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
                plt.show(block = False)

            in_text = raw_input("Input: \n" +
                                "\"DONE\" to stop\n" +
                                "\"D#\" for a vertex index to delete. \n" +
                                "\"C#,#\" to move the first vertex to the location of the second.\n" +
                                "\"UPDATE\" to refresh the screen.")

            if in_text == "UPDATE":
                refresh = True
                continue
            if in_text == "DONE":
                break
            try:
                if in_text.startswith("C"):
                    first_vert_txt, second_vert_txt = in_text.split(',')
                    first_vert = int(first_vert_txt[1:])
                    second_vert = int(second_vert_txt)
                    combine_list.append((first_vert, second_vert))
                    twod_vertices[first_vert] = twod_vertices[second_vert]
                elif in_text.startswith("D"):
                    vert_text = in_text[1:]
                    del_vert = int(vert_text)
                    del_vert_list.append(del_vert)
                    del twod_vertices[del_vert]
            except:
                print("Bad input. Try again.")
                continue
    # FOR:
    # basin_pt = [29.439518, 106.485528]
    # tibet_pt = [34.7327312, 100.8766929]
    del_vert_list = [45,44,63,62,74,78,36,34,35,70,0]
    combine_list = [(64,1)]
    # FOR:
    # basin_pt = [29.5, 106.5]
    # tibet_pt = [34,101]
    # del_vert_list = [23, 22, 55, 80, 81, 67, 51, 50, 66]
    # combine_list = [(79, 0)]
    d['del_vert_list'] = del_vert_list
    d['combine_list'] = combine_list

def process_fixes(d):
    for v1, v2 in d['combine_list']:
        d['twod_vertices'][v1] = d['twod_vertices'][v2]

    # Sort in reverse order before deleting so we don't screw up the
    # indexing.
    for vert in d['del_vert_list']:
        del d['twod_vertices'][vert]

    d['twod_segments'] = [seg for seg in d['twod_segments']
                          if seg[0] in d['twod_vertices']
                             and seg[1] in d['twod_vertices']]

