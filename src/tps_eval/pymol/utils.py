def show_organic_and_metals(cmd):
    cmd.show("sticks", "structure and organic")
    cmd.color("atomic", "structure and organic")
    cmd.color("gray", "structure and organic and elem C")
    cmd.show("spheres", "structure and metals")
