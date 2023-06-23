from math import *
import matplotlib.pyplot as plt
import time
import random as rand
from math import log10

#Hellooooo

## Classes
class Point:
    def __init__(self,x,y):
        self.x = x
        self.y = y
    def equals(self, other):
        return (self.x == other.x) and (self.y == other.y)
    def __str__(self):
        return "Point("+str(self.x)+", "+str(self.y)+")"

class Vecteur:
    def __init__(self,dx,dy):
        self.dx = dx
        self.dy = dy

    def produitScalaire(self, other):
        return self.dx*other.dx + self.dy*other.dy

    def normal(self):
        return Vecteur(-self.dy, self.dx)

    def norme(self):
        return sqrt( self.produitScalaire(self) )

    def multiplie(self, valeur):
        self.dx *= valeur
        self.dy *= valeur

    def devientUnitaire(self):
        maNorme = self.norme()
        self.dx /= maNorme
        self.dy /= maNorme

    def __str__(self):
        return ">"+str(self.dx)+", "+str(self.dy)+">"

class Droite:
    def __init__(self,pt1,pt2):
        if (pt1.equals(pt2)):
            raise Exception("BUG: on ne peut pas definir une doite avec 2x le meme point")
        self.v = Vecteur((pt2.x - pt1.x),(pt2.y - pt1.y))
    #un vecteur directeur de la droite
        self.a = self.v.dy
        self.b = -self.v.dx
        self.c = (-self.a*pt1.x -self.b*pt1.y)
    #coefficients dans l'équation de droite ax + by + c = 0
        self.pt1 = pt1
        self.pt2 = pt2
    #renvoie un point appartenant à la droite

    def __str__(self):
        return "Droite{pt1="+str(self.pt1)+", pt2="+str(self.pt2)+"}"

    def intersection(d1,d2): #Renvoie l'intersection
        if d2.a*d1.b - d2.b*d1.a == 0 : #si a2*b1-b2*a1 = 0, les droites sont parallèles
            return None
        else:
            x = (d2.b * d1.c - d1.b * d2.c) / (d2.a * d1.b - d2.b * d1.a)
            y = (d2.a * d1.c - d1.a * d2.c) / (d1.a * d2.b - d1.b * d2.a)
            return Point(x,y)

class Site :
    def __init__(self,pt,id):
        self.pt = pt
        self.id = id
        self.adja = {}
        #self.d_virt = d_virt #None si c'est un vrai site, la droite de la fenêtre correspondant au site virtuel sinon
    def distance(self, other):
        return distance_site(self, other)
    def __str__(self):
        return "Site#"+str(self.id)+"{"+str(self.pt)+"}"

class Edge:
    def __init__(self,pt1,pt2):#,site1,site2):
        if (pt1.equals(pt2)):
            raise Exception("BUG: on ne peut pas definir une doite avec 2x le meme point")
        self.pt1 = pt1
        self.pt2 = pt2
    def droite(self):
        return Droite(self.pt1, self.pt2)
    def milieu(self):
        return Point((self.pt1.x+self.pt2.x)/2, (self.pt1.y+self.pt2.y)/2)
    def __str__(self):
        return "["+str(self.pt1)+", "+str(self.pt2)+"]"


class Fenetre:#rajouter les points virtuels
    def __init__(self,x1,x2,y1,y2):
        self.xmin = min(x1, x2)
        self.xmax = max(x1, x2)
        self.ymin = min(y1, y2)
        self.ymax = max(y1, y2)
        ptGB = Point(self.xmin,self.ymin) #en bas/gauche
        ptGH = Point(self.xmin,self.ymax) #en haut/gauche
        ptDH = Point(self.xmax,self.ymax) #en haut/droite
        ptDB = Point(self.xmax,self.ymin) #en bas/droite
        eG = Edge(ptGB,ptGH) #gauche
        eH = Edge(ptGH,ptDH) #haut
        eD = Edge(ptDH,ptDB) #droite
        eB = Edge(ptDB,ptGB) #bas
        self.bords = [ eG, eB, eD, eH ]

class DiagVoronoi:
    def __init__(self,liste_sites,fenetre): #les 4 sites virtuels sont les 4 premiers sites de liste_sites
        self.lsites = liste_sites
        self.fenetre = fenetre
        self.debugGraphId = 0

def dumpSite(site):
    msg = str(site)
    for siteAdja in site.adja.values():
        msg += "\n\t-> "
        msg += str(siteAdja)
    return msg

## Fonctions géométriques

def distance_pt(pt1,pt2):
    return sqrt((pt2.x - pt1.x)**2 + (pt2.y - pt1.y)**2)

def distance_site(siteI,siteJ):
    return distance_pt(siteI.pt,siteJ.pt)

def medIJ(siteI,siteJ):
    milieu = Point( (siteI.pt.x + siteJ.pt.x)/2, (siteI.pt.y + siteJ.pt.y)/2 ) #pt au milieu de siteI et siteJ. €g
    d = Droite(siteI.pt,siteJ.pt) #droite entre I et J, perpendiculaire à la médiatrice de I et J
    v = d.v
    n = v.normal() #vecteur normal à dIJ, directeur de g
    pt2 = Point( (milieu.x + n.dx), (milieu.y + n.dy)) # autre point de g
    return Droite(milieu,pt2) #g


def interDrSeg(d,s): #s est une arête (edge, segment) et d une droite
    Ds = Droite(s.pt1,s.pt2)
    pt_inter = d.intersection(Ds)

    ##print("interDrSeg:")
    ##print("\t* droite = "+str(d))
    ##print("\t* segmnet = "+str(s))
    ##print("\t* pt_inter = "+str(pt_inter))

    precision = 0.0001

    if pt_inter == None :
        return None

    if pt_inter.x > s.pt1.x+precision and pt_inter.x > s.pt2.x+precision :
        return None

    if pt_inter.x < s.pt1.x-precision and pt_inter.x < s.pt2.x-precision :
        return None

    if pt_inter.y > s.pt1.y+precision and pt_inter.y > s.pt2.y+precision :
        return None

    if pt_inter.y < s.pt1.y-precision and pt_inter.y < s.pt2.y-precision :
        return None

    #tests pour voir si pt_inter est dans le segment

    return pt_inter

##Green Sibson

##########################################################################
##########################################################################
##########################################################################

def nearest_site(siteN,siteDep):#siteDep est un site du diagramme déjà existant
    siteMin = siteDep
    dstMin  = distance_site(siteN, siteMin)
    while True:
        for siteAdja in siteDep.adja.keys():
            dstAdja = siteN.distance(siteAdja)
            if dstAdja < dstMin :
                dstMin = dstAdja
                siteMin = siteAdja

        if siteMin == siteDep:
            return siteDep

        # OPTIM : on evite un appel recursif de fonction (+ couteux que juste un saut vers le debut de la boucle, ne fonctionne pas si trop de niveau de recursion)
        siteDep = siteMin

##########################################################################
##########################################################################
##########################################################################

##########################################################################
# renvoie si pt1 et pt2 sont du meme cote de la droite                   #
##########################################################################
def memeCote(droite, pt1, pt2):
    v1 = Vecteur(droite.pt1.x-pt1.x, droite.pt1.y-pt1.y)
    v2 = Vecteur(droite.pt1.x-pt2.x, droite.pt1.y-pt2.y)
    vn = droite.v.normal()
    # astuce : si les 2 produits scalaire sont du meme signe, alors les 2 pts sont du meme cote de la droite
    #          meme signe <=> le produit et positif
    return (vn.produitScalaire(v1) * vn.produitScalaire(v2)) > 0

##########################################################################
# ajuste l'edge commun aux sites I et J suite a l'apparition du site K   #
##########################################################################
def modifie_edge(media, siteI, siteJ, new_pt, diag): #new_pt est le point d'intersection de l'arête N-J avec I-J. L'arête I-J doit être mutilée en partie
    ##print("---------------------------------------")
    ##print("ENTER : modifie_edge(")
    ##print("\tmedia    : "+str(media))
    ##print("\tsiteI    : "+str(siteI))
    ##print("\tsiteJ    : "+str(siteJ))
    ##print("\tnew_pt   : "+str(new_pt))
    ##print(")")
    ##print("---------------------------------------")
    edge = siteJ.adja[siteI]
    # astuce : on fait une copie, car on veut garder les anciens coordonnees visible dans le sens
    #edge = Edge(edge.pt1, edge.pt2)
    #siteJ.adja[siteI] = edge
    # optim : plutot que de recreer un Edge, on modifie l'ancien edge deja existant
    if memeCote(media, siteJ.pt, edge.pt1):
        edge.pt2 = new_pt
    else:
        edge.pt1 = new_pt
    debugGraph(diag)
    #siteTronque.adja[siteK] = eIJ
    #siteK.adja[siteTronque] = eIJ

###########################################################################################
# description:                                                                            #
#   cherche la 1ier arrete qui se fait couper par la mediatrice separant les sites K et N #
# entrees:                                                                                #
#   siteN : site qu'on ajoute                                                             #
#   siteK : site impacte par siteN                                                        #
#   siteLePlusProche : le + proche des site impacte par siteN                             #
#   diag  : le diagramme                                                                  #
# sortie:                                                                                 #
#   id du 1ier site adjacent trouve avec dont l'edge se fait couper                       #
###########################################################################################
def first_edge(siteN, siteK, diag):

    ##print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    ##print("ENTER: first_edge("+str(siteN)+", "+str(siteK)+")")
    ##print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

    if siteN in siteK.adja:
        return # on a deja travailler sur ce cas. ce test permet d'evitre de boucler indefiniment

    sitesTronques = {} # clef = site, valeur = edge
    sitesSupprimes = []
    sitesGardes = 0

    media = medIJ(siteN, siteK)
    for siteI in siteK.adja.keys():
        eIK = siteK.adja[siteI]
        if memeCote(media, eIK.pt1, siteK.pt):
            if memeCote(media, eIK.pt2, siteK.pt):
                # eIK est entierement du meme cote de la mediatrice que le centre sur siteK
                # eIK n'est ni tronquee, ni supprimee
                sitesGardes += 1
                continue
            # eIK est a tronquer
            sitesTronques[siteI] = eIK
            continue
        if memeCote(media, eIK.pt2, siteK.pt):
            # eIK est a tronquer
            sitesTronques[siteI] = eIK
            continue
        # eIK est entierement de l'autre cote de la mediatrice, les sites I et K se deconnectent
        sitesSupprimes.append(siteI)

    # on deconnecte asap tous les sites qui sont impactes par l'insert du site N
    for siteI in sitesSupprimes:
        del siteK.adja[siteI]
        debugGraph(diag)

    # si pas de bug, la mediatrice coupe 2 et uniquement 2 edge
    # si on a un seul tronque (ou 0), c'est que l'autre tronque est un des bords. TODO: trouver le bord pour calcul l'arrete
    if len(sitesTronques) > 2:
        raise Exception("WTF ? on a "+str(len(sitesTronques))+" edges a tronquer en ajoutant "+str(siteN)+" dans "+str(siteK)+". Pour infos, on avait "+str(sitesGardes)+" sites gardes et "+str(len(sitesSupprimes))+" sites supprimes")

    pts = []
    if len(sitesTronques) == 0:
        # la mediatrice coupe 2 bords et seulement 2 bords (on ne traite pas le cas ou on arrive exactement dans 1 angle :/)
        #print("\taucun edge tronque, on cherche 2 bords")
        #print("\tpour rappel, la mediatrice est "+str(media))
        for bord in diag.fenetre.bords:
            #print("\t\ton test avec le bord "+str(bord))
            #print("\t\tpour rappel, la mediatrice est "+str(media))
            pt = interDrSeg(media, bord)
            #print("\t\ttouche la mediatrice en "+str(pt))
            if pt != None:
                pts.append(pt)
                if len(pts)==2:
                    break
        #print("\ton a trouve "+str(len(pts))+" intersection bord/mediatrice sur les 2 attendus")
    else:
        # pas 0, donc soit 1 ou 2
        # dans les 2 cas, la ou les arretes tronquees sont a ajuster
        for siteTronque in sitesTronques.keys():
            pt = interDrSeg(media, siteK.adja[siteTronque])
            pts.append( pt )
            modifie_edge(media, siteTronque, siteK, pt, diag)

        # et si seulement 1, on essaye de trouver l'autre point sur l'un des bords
        if len(sitesTronques) == 1:
            #print("\ton a un seul edge tronque, donc il faut trouver le bord qu'on coupe")
            # on touche un bord, certe, mais ce doit etre un bord que touche aussi l'autre site
            for bord in diag.fenetre.bords:
                #print("\t\ton essaye avec le bord "+str(bord))
                ptBord = interDrSeg(media, bord)
                #print("\t\tqui coupe la mediatrice en "+str(ptBord))
                if ptBord == None:
                    continue
                # ok, on a une intersection, mais elle est peut etre du mauvais cote
                # pour cela, on test tous les milieux de tous les edges restant du site
                # qui doivent tous etre du bon cote (la zone d'influence est convexe)
                # donc le point du bord et le point du site doivent tjs etre du
                # meme code que les droites formees par les edges actuellement connu
                pasLeBonBord = False
                for edge in siteK.adja.values():
                    droite = edge.droite()
                    if not memeCote(droite, siteK.pt, ptBord):
                        pasLeBonBord = True
                        break
                if pasLeBonBord:
                    continue
                #print("!!!! trouve !!!!!")
                break

            pts.append( ptBord )

    # on rajoute la nouvelle arrete complete entre les site N et K
    #print("on avait "+str(len(sitesTronques))+" edges intersites tronque:")
    #for pt in pts:
        #print("\t"+str(pt))
    siteK.adja[siteN] = Edge(pts[0], pts[1])
    siteN.adja[siteK] = Edge(pts[0], pts[1])
    debugGraph(diag)

    # on fait un appel recursif sur tous les sites supprime. Ils ne peuvent pas nous rappeler car on a deja deconnecter siteK de tous ces sites supprimes un peu + haut
    for siteSupprime in sitesSupprimes:
        first_edge(siteN, siteSupprime, diag)

    for siteTronque in sitesTronques.keys():
        first_edge(siteN, siteTronque, diag)


###########################################################################
# description                                                             #
#   ajoute un nouveau site encore inconnue a cote du seul site existant   #
#   qui l'inclue dans sa zone d'influence                                 #
# entrees                                                                 #
#   siteN : le nouveau site                                               #
#   siteK : le site accueillant                                           #
#   diag  : le diagramme                                                  #
# sorties                                                                 #
#   aucune                                                                #
###########################################################################
def ajout_nouveau_site_a_son_plus_proche(siteN, siteK, diag):
    first_edge(siteN, siteK, diag)

def GreenSibson__V1(diag):
    site_prec = diag.lsites[0]
    for siteN in diag.lsites[1:]:
        ##print("==============================================================")
        ##print(f"GreenSibson__V1, prochain site a ajouter: {siteN}")
        ##print("==============================================================")
        siteK = nearest_site(siteN,site_prec)
        ##print(f"\tSite le + proche: {siteK}")
        ##print(dumpSite(siteK))
        ajout_nouveau_site_a_son_plus_proche(siteN, siteK, diag)
        site_prec = siteN
    return diag

GreenSibson = GreenSibson__V1

def graph(sites, fenetre, filename):
    plt.figure()

    # on dessine tous les sites
    x_sites = []
    y_sites = []
    for site in sites:
        x_sites.append(site.pt.x)
        y_sites.append(site.pt.y)
    plt.plot(x_sites,y_sites,color='black',ls='None',marker='+')

    # on dessine les contours des sites
    for site in sites:
        for edge in site.adja.values() :
            x_aretes = []
            y_aretes = []
            x_aretes.append(edge.pt1.x)
            x_aretes.append(edge.pt2.x)
            y_aretes.append(edge.pt1.y)
            y_aretes.append(edge.pt2.y)
            plt.plot(x_aretes,y_aretes,color='black',marker='None',ls='-')

    # on dessine aussi les bords
    x_aretes = []
    y_aretes = []
    for edge in fenetre.bords :
        x_aretes.append(edge.pt1.x)
        y_aretes.append(edge.pt1.y)
    edge = fenetre.bords[0]
    x_aretes.append(edge.pt1.x)
    y_aretes.append(edge.pt1.y)
    plt.plot(x_aretes,y_aretes,color='red',marker='None',ls='-')

    plt.savefig(filename)
    plt.show()

def debugGraph(diag):
    return
    #graph(diag.lsites, diag.fenetre, "debug_"+str(diag.debugGraphId)+".png")
    #diag.debugGraphId += 1

## Tests

def generateur_aleat(n, f):
    #rand.seed(13)
    #ls = [f.siteG,f.siteH,f.siteD,f.siteB]
    ls = []

    epsx = (f.xmax - f.xmin)/100
    epsy = (f.ymax - f.ymin)/100

    for k in range(n):
        sx=rand.uniform(f.xmin + epsx, f.xmax - epsx)
        sy=rand.uniform(f.ymin + epsy, f.ymax - epsy)
        #insert(ls, Site( Point(sx,sy), k+4, None ) )
        ls.append(Site( Point(sx,sy), k ))

    return ls


fenetre = Fenetre(0, 100, 0, 100)
sites = generateur_aleat(10, fenetre)

##print("Test calculs intersections droite/segment -------------------------")
media = Droite(Point(8.918023245394973, 25.548989685562198), Point(13.884067246345293, 38.55628715578214))
media.v.devientUnitaire()
nbInter = 0
#print("droite de test: "+str(media))
#print("vecteur directeur de cette droite: "+str(media.v))
for bord in fenetre.bords:
    pt = interDrSeg(media, bord)
    #print("l'intersection entre "+str(media)+" avec "+str(bord)+" est "+str(pt))
    if pt != None:
        nbInter += 1
        v = Vecteur(pt.x-media.pt1.x, pt.y-media.pt1.y)
        v.devientUnitaire()
        #print("vecteur directeur de l'intersection: "+str(v))
if nbInter != 2:
    raise Exception("Bug dans le calcul d'intersection droite/segment detecte. On devrait avoir 2 points, on en a trouve "+str(nbInter))
#print("------------------------- Test calculs intersections droite/segment")

#print("Les sites ---------------------------------")
#for site in sites:
    #print(site)
#print("--------------------------------- Les sites")

##print("Test distances au 4 bords -----------------")
##print(distance_site(fenetre.siteG, sites[4]))
##print(distance_site(sites[4], fenetre.siteG))
##print(distance_site(fenetre.siteH, sites[4]))
##print(distance_site(sites[4], fenetre.siteH))
##print(distance_site(fenetre.siteD, sites[4]))
##print(distance_site(sites[4], fenetre.siteD))
##print(distance_site(fenetre.siteB, sites[4]))
##print(distance_site(sites[4], fenetre.siteB))
##print("----------------- Test distances au 4 bords")

#print("Test distances 2 pts ----------------------")
#print(distance_site(sites[4], sites[5]))
#print(distance_site(sites[5], sites[4]))
#print("---------------------- Test distances 2 pts")

#print("Test produit scalaire ---------------------")
media = medIJ(sites[4], sites[5])
v = Vecteur( sites[4].pt.x-sites[5].pt.x, sites[4].pt.y-sites[5].pt.y )
#print(v.produitScalaire(v.normal()))
#print(media.v.produitScalaire(Vecteur( sites[4].pt.x-sites[5].pt.x, sites[4].pt.y-sites[5].pt.y )))
#print("--------------------- Test produit scalaire")

tab_pts = [i for i in range(1,500)]
stats = []
nbPts = 1
#print("nbPts = "+str(nbPts))
while nbPts<500:
    #print("nbPts = "+str(nbPts))
    retry = 20
    stats_current=[]
    while retry>0:
        retry -= 1
        sites = generateur_aleat(nbPts, fenetre)
        diagtest1 = DiagVoronoi(sites,fenetre)
        ts = time.time()
        GreenSibson(diagtest1)
        ts = time.time()-ts
        stats_current.append(ts)
        #graph(sites, fenetre, "zePlot_"+str(nbPts)+".png")

    nbPts += 1
    stats_current.sort()
    mediane = (stats_current[9] + stats_current[10])/2
    stats.append(mediane)
    #plt.show()


plt.figure()
plt.plot(tab_pts,stats)
plt.show()

