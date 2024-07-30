def cad():
    import pandas as pd
    from shutil import copyfile

    
    df = pd.read_excel("cad.xlsx")
    
    #n = len(df.index) - 1
    copyfile("cylshell.txt", "cad.txt")
    fi = open("cad.txt","a+")
    M = 0
    xcm = 0
    ycm = 0
    zcm = 0
    for index, row in df.iterrows():

        x = row['xc']
        y = row['yc']
        z = row['zc']
        ty = row['type']
        m = row['m']
        M = M+m
        xcm = xcm + x*m
        ycm = ycm + y*m
        zcm = zcm + z*m
        
        if ty == "cube":
            a = row['a']
            b = row['b']
            c = row['c']
            str = "translate([{},{},{}]) cube([{},{},{}]);\n".format(x,y,z,a,b,c)
        elif ty == "cylinder":
            r = row['r']
            h = row['h']
            str = "translate([{},{},{}]) cylinder(r={},h={},$fn=20);\n".format(x,y,z,r,h)
        elif ty == "sphere":
            r = row['r']
            str = "translate([{},{},{}]) sphere({},$fn=20);\n".format(x,y,z,r)
        elif ty == "cylshell":
            r = row['r']
            h = row['h']
            t = row['t']
            str = "translate([{},{},{}]) cylshell({},{},{});\n".format(x,y,z,r,h,t)
        elif ty == "sphershell":
            r = row['r']
            t = row['t']
            str = "translate([{},{},{}]) sphershell({},{});\n".format(x,y,z,r,t)
        fi.write(str)
    xcm = xcm / M
    ycm = ycm / M
    zcm = zcm / M
    fi.close()
    return xcm,ycm,zcm