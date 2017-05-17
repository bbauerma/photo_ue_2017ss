import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

if __name__ == "__main__":

    np.set_printoptions(precision=6)
    np.set_printoptions(suppress=True)

    plt.style.use("ggplot")

    #=======================================================================
    # test set used for controlling the results of the calculation of the
    # fundamentalmatrix; points can be found in KRAUSS - PHOTOGRAMMETRIE BAND 1
    # CHAPTER 4.3.1 - Page 219
    # ======================================================================
    #p1 = np.array([
        #[1,93.176,5.890],
        #[2,-27.403,6.672],
        #[3,83.951,107.422],
        #[4,-11.659,101.544],
        #[5,110.326,-97.800],
        #[6,-12.653,-87.645],
        #[7,37.872,40.969],
        #[8,41.503,-37.085]
        #]
    #)

    #p2 = np.array([
        #[1,6.072,5.176],
        #[2,-112.842,1.121],
        #[3,-4.872,105.029],
        #[4,-99.298,95.206],
        #[5,34.333,-99.522],
        #[6,-96.127,-93.761],
        #[7,-48.306,37.862],
        #[8,-42.191,-40.138]
        #]
    #)

    #camera_const = 152.67 #mm
    #=======================================================================
    #EXERCISE VALUES
    p1 = np.array([
        [54, 1505.791, -934.224],
        [81, 794.0793, -943.624],
        [85, 1429.979, -1219.010],
        [129, 831.493, -1154.156],
        [159, 1722.855, -1221.836],
        [193, 889.724, -1135.457],
        [215, 1778.163, -1298.671],
        [237, 1536.589, -901.604],
        [361, 1049.607, -514.139],
        [407, 1944.281, -1256.023],
        [468, 1639.276, -1051.416],
        [510, 1867.398, -716.761]
    ])


    p2 = np.array([
        [54, 1893.071, -718.891],
        [81, 1217.832, -1371.529],
        [85, 1708.045, -1159.236],
        [129, 1159.775, -1541.060],
        [159, 1960.698, -894.555],
        [193, 1105.454, -1529.853],
        [215, 1912.039, -981.543],
        [237, 1594.205, -829.155],
        [361, 1024.032, -1060.022],
        [407, 1730.554, -915.923],
        [468, 1227.926, -1121.613],
        [510, 1505.277, -458.463]
    ])

    camera_const = float(3400)  #px
    #===================================================================

    img_width = 3000  #px
    img_height = 2000 #px
    origin = [0.5, -0.5] #px

    #PRETTY PRINTING OF THE TIEPOINTS; NOTHING ELSE...
    print "%i TIEPOINTS: " % np.shape(p1)[0]
    print "------------------------------------------------------------"
    for i, x in enumerate(p1):

        if i == 0:
             print "{:>5}".format("ID") + " " + "{:>11}".format("P1_X") + " " + "{:>11}".format("P1_Y") + " " + "{:>14}".format("P2_X") + " " + "{:>11}".format("P2_Y")

        p_id = p1[i,0]

        p1_x = p1[i,1]
        p1_y = p1[i,2]

        p2_x = p2[i,1]
        p2_y = p2[i,2]



        print "{:>5}".format(str(p_id)) + " - " + "{:10.4f}".format(p1_x) + " " + "{:10.4f}".format(p1_y) + "  --> " + "{:10.4f}".format(p2_x) + " " + "{:10.4f}".format(p2_y)

    print "------------------------------------------------------------"
    print ""
    # ==================================================================



    principal_point = [(img_width - 1) / 2., -(img_height - 1) / 2.]
    calib_matrix = np.array([
                    [1,0,-principal_point[0]],
                    [0,1,-principal_point[1]],
                    [0,0,-camera_const]
                    ]
    )


    ## =================================================================
    ## CALCULATION OF F MATRIX
    ## =================================================================
    mat_a = np.zeros((np.shape(p1)[0], 9))

    #A contains kron(P2,P1)^T in each row (kronecker product);
    #P1 and P2 must be in shape (x,y,1)^T --> homogeneous coords
    for i,x in enumerate(p1):

        p1_vec = np.array([[p1[i,1]], [p1[i,2]], [1]])
        p2_vec = np.array([[p2[i,1]], [p2[i,2]], [1]])

        mat_a[i,:] = np.transpose(np.kron(p2_vec, p1_vec))

    mat_at = np.transpose(mat_a)
    mat_ata = np.dot(mat_at, mat_a)

    #SVD - SINGLE VALUE DECOMPOSITION
    #fundamental matrix f is the eigenvector corresponding to the smallest
    #eigenvalue in s; as the eigenvalues in s are sorted DESC the smallest
    #is s[-1]; the resultung fundamentalmatrix f in general won't have
    #det(F) == 0

    U1,s1,V1 = np.linalg.svd(mat_ata)

    f = V1[-1]
    mat_F = np.transpose(np.reshape(f, (3,3)))

    #normalize the last element to 1; not necessary but done in KRAUSS
    mat_F = mat_F / mat_F[2,2]

    print "\n\nFUNDAMENTALMATRIX (det!=0) = "
    print "------------------------------------------------------------"
    print mat_F
    print "\nDET(F) = %.6f" % np.linalg.det(mat_F)

    print "\nWIDERSPRUCH P1 * F * P2 != 0"
    for i,x in enumerate(p1):

        p1_vec = np.array([[p1[i,1]], [p1[i,2]], [1]])
        p2_vec = np.array([[p2[i,1]], [p2[i,2]], [1]])

        print np.dot(np.dot(np.transpose(p1_vec), mat_F),p2_vec)

    print "------------------------------------------------------------"
    print ""

    #because in general F won't have det(F) = 0 we need to SVD our original
    #fundamentalmatrix and change the smallest eigenvalue to 0; by
    #calculating U*S*V we get F' which should have a det(F') of 0;
    U2,s2,V2 = np.linalg.svd(mat_F)
    S2 = np.zeros(np.shape(mat_F))

    S2 = np.diag(s2)

    S3 = np.copy(S2)
    S3[2,2] = 0

    #F1 is wrong at the moment....everything else works fine
    mat_F1 = np.dot(np.dot(U2, S3), V2)

    print "\n\nFUNDAMENTALMATRIX F' (det == 0) = "
    print "------------------------------------------------------------"
    print mat_F1
    print "\nDET(F') = %.6f" % np.linalg.det(mat_F1)

    print "\nWIDERSPRUCH P1 * F' * P2 != 0 "
    for i,x in enumerate(p1):

        p1_vec = np.array([[p1[i,1]], [p1[i,2]], [1]])
        p2_vec = np.array([[p2[i,1]], [p2[i,2]], [1]])

        print np.dot(np.dot(np.transpose(p1_vec), mat_F1),p2_vec)

    print "------------------------------------------------------------"

    # ==================================================================
    #CALCULATION OF ROTATION MATRIX R1 AND R2
    #V = A*X-l
    #v = mat_a * x - l_obs
    # ==================================================================

    #dim(mat_a) = [nr of tiepoints; 5] #dk1, dk2, dph1, dph2, dom2
    #mat_a = np.zeros((np.shape(p1)[0],5))

    #mat_a[:,0] = -p1[:,1]
    #mat_a[:,1] = p2[:,1]
    #mat_a[:,2] = (p1[:,1]*p1[:,2]) / camera_const
    #mat_a[:,3] = -(p2[:,1]*p2[:,2]) / camera_const
    #mat_a[:,4] = camera_const + (p2[:,2]**2 / camera_const)

    #l_obs = np.subtract(p1[:,2] , p2[:,2])

    #mat_at = np.transpose(mat_a)

    #v_rad = np.dot(np.dot(np.linalg.inv(np.dot(mat_at, mat_a)), mat_at), l_obs)
    #v_gon = v_rad*200/np.pi
    # ==================================================================


    # ==================================================================
    # PLOTTING THE IMAGES AND ANYTHING ELSE
    # ==================================================================
    #min_x_p1 = np.min(p1[:,1])
    #max_x_p1 = np.max(p1[:,1])
    #min_y_p1 = np.min(p1[:,2])
    #max_y_p1 = np.max(p1[:,2])

    #min_x_p2 = np.min(p2[:,1])
    #max_x_p2 = np.max(p2[:,1])
    #min_y_p2 = np.min(p2[:,2])
    #max_y_p2 = np.max(p2[:,2])

    #img_bbox_x = [0, 3000, 3000, 0,0]
    #img_bbox_y = [0, 0, -2000, -2000,0]

    ##min_x = min([min_x_p1, min_x_p2]) - (abs(max_x) - abs(min_x)) / 10
    ##max_x = min([max_x_p1, max_x_p2]) + (abs(max_x) - abs(min_x)) / 10
    ##min_y = min([min_y_p1, min_y_p2]) - (abs(max_y) - abs(min_y)) / 10
    ##max_y = min([max_y_p1, max_y_p2]) + (abs(max_y) - abs(min_y)) / 10

    #min_x = -100
    #max_x = 3100
    #min_y = -2100
    #max_y = 100

    ## Two subplots, unpack the axes array immediately
    #f, (ax1, ax2) = plt.subplots(1, 2, figsize=(30,10), sharey="all")

    ##plot image frame
    #ax1.add_patch(
        #patches.Rectangle(
            #(0, 0), # (x,y)
            #3000,   # width
            #-2000, # height
            #facecolor=(0,0,0,0.1),
            #edgecolor="black",
            #linewidth=1,
        #)
    #)

    #ax2.add_patch(
        #patches.Rectangle(
            #(0, 0), # (x,y)
            #3000,   # width
            #-2000, # height
            #facecolor=(0,0,0,0.1),
            #edgecolor="black",
            #linewidth=1,
        #)
    #)

    ##IMAGE P1
    #ax1.scatter(p1[:,1], p1[:,2], s=100, facecolors='red', alpha=.5, edgecolors='r', linewidth=1.5)
    #ax1.scatter(principal_point[0], principal_point[1], marker="x", s=100, facecolor="black", edgecolors="black", linewidth=2)

    #for i, txt in enumerate(p1[:,0]):
        #ax1.annotate(int(txt), (p1[i,1]+30,p1[i,2]-30), fontsize=15, fontweight="bold")

    ##IMAGE P2
    #ax2.scatter(p2[:,1],p2[:,2], s=100, facecolors='blue', alpha=.5, edgecolors='b', linewidth=1.5)
    #ax2.scatter(principal_point[0], principal_point[1], marker="x", s=100, facecolor="black", edgecolors="black", linewidth=2)

    #for i, txt in enumerate(p2[:,0]):
        #ax2.annotate(int(txt), (p2[i,1]+30,p2[i,2]-30), fontsize=15, fontweight="bold")

    ##SET AXIS RANGE
    #ax1.set_xlim([min_x, max_x])
    #ax1.set_ylim([min_y, max_y])
    #ax2.set_xlim([min_x, max_x])
    #ax2.set_ylim([min_y, max_y])

    ##PLOT OPTIONS
    #plt.tight_layout()
    ##plt.show()
    #plt.savefig("./photo_2_points.png")
