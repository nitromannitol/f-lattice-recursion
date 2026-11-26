# Builds the tiles for a word one letter at a time.

_recphase(w) = count(==('1'), w) % 2
_recall(w,c) = all(==(c),w)
function _recword(w::AbstractString)
    all(c -> c in ('1','2','3'),w) || throw(ArgumentError("word must contain only 1, 2, 3"))
    String(w)
end
function _recfirst(w)
    isempty(w) && return ""
    w[1]=='1' && return ""
    c=w[1]; i=1
    while i<=length(w) && w[i]==c; i+=1; end
    w[1:i-1]
end
function _recmeld(tiles, order; context="")
    out=copy(tiles[first(order)])
    for j in Iterators.drop(order,1)
        out=touch(out,tiles[j]; context="$context slot$(j-1)")
    end
    trim(out)
end
function _recplaced(tiles,vs,slopes)
    [place(tiles[i],vs[i],slopes[i]) for i in 1:length(tiles)]
end

# Find the first prefix of w at which frac appears as the odd (or even) fraction.
function _recbirth(frac,w,odd)
    slot=odd ? 1 : 2
    for k in 0:length(w)
        a=w[1:k]
        fs=getFracs(a)
        if fs[slot]==frac
            return a,fs
        end
    end
    error("$frac does not appear along the word $w")
end
function translate_alts(mats,w)
    palt,qalt=mats[3],mats[4]
    p,q=getFracs(w)[3:4]
    if q != 1//1
        birth,fs=_recbirth(q,w,false)
        if _recphase(birth)==0
            O=vf(fs[4])[1]-vf(fs[3])[2]
            if imag(O)>=0
                _recall(birth,'3') && (O+=pt(0,1))
                qalt=place(qalt,-O,pt(0,0))
            end
        end
    end
    if p != 0//1
        birth,fs=_recbirth(p,w,true)
        if _recphase(birth)==1
            O=vf(fs[4])[1]-vf(fs[3])[2]
            imag(O)>=0 && (palt=place(palt,-O,pt(0,0)))
        end
    end
    palt,qalt
end
function get_std_child_unc(mats,w;getOdd=true)
    p,q=mats[1:2]
    vs,aa=getOdd ? getOddTileFull(w) : getEvenTileFull(w)
    tiles=getOdd ? (q,q,q,q,p,p) : (q,q,p,p,p,p)
    order=getOdd ? (1,2,5,6,3,4) : (1,3,4,2,5,6)
    _recmeld(_recplaced(tiles,vs,aa),order; context="$w standard unc odd=$getOdd")
end
function get_alt_child_unc(mats,w;getOdd=true)
    p,q=mats[1:2]
    pa,qa=translate_alts(mats,w)
    vs,aa=getOdd ? getAltOddTileFull(w) : getAltEvenTileFull(w)
    tiles=getOdd ? (q,q,p,p,p,p) : (q,q,q,q,p,p)
    placed=_recplaced(tiles,vs,aa)
    v,a=get_alt_full_offset(getFracs(w)[getOdd ? 1 : 2],w)
    push!(placed,place(getOdd ? qa : pa,v,a))
    order=getOdd ? (1,7,3,4,2,5,6) : (1,7,2,5,3,4,6)
    _recmeld(placed,order;context="$w alternate unc odd=$getOdd")
end
function _recalign(mats,helpers,w,standardOdd)
    gl,gr=helpers
    target=ur(mats[standardOdd ? 2 : 1])
    if w[end]!='3' && _recphase(w)==1
        gr=place(gr,target-ur(gr),pt(0,0))
    else
        gl=place(gl,target-ur(gl),pt(0,0))
    end
    gl,gr
end
function get_std_child_c(mats,helpers,w;getOdd=false)
    p,q=mats[1:2]; gl,gr=_recalign(mats,helpers,w,getOdd)
    vs,aa=getOdd ? getOddTileFull(w) : getEvenTileFull(w)
    tiles=getOdd ? (_recphase(w)==0 ? (q,gr,gl,q,p,p) : (gr,q,q,gl,p,p)) : (q,q,p,gl,gr,p)
    order=getOdd ? (1,2,5,6,3,4) : (1,3,4,2,5,6)
    _recmeld(_recplaced(tiles,vs,aa),order;context="$w standard corrected odd=$getOdd")
end
function get_alt_child_c(mats,helpers,w;getOdd=true)
    p,q=mats[1:2]; gl,gr=_recalign(mats,helpers,w,!getOdd)
    vs,aa=getOdd ? getAltOddTileFull(w) : getAltEvenTileFull(w)
    tiles=getOdd ? (q,q,p,gl,gr,p) : (gr,q,q,gl,p,p)
    _recmeld(_recplaced(tiles,vs,aa),(1,3,4,5,6,2);context="$w alternate corrected odd=$getOdd")
end

function get_l_correct(mats,alts,ancestor,w;getOdd=true,isFlip=false,isLong=false)
    k=length(w)-length(ancestor)
    startswith(w,ancestor) || throw(ArgumentError("ancestor is not a prefix of the word"))
    all(==(getOdd ? '2' : '3'),w[length(ancestor)+1:end]) ||
        throw(ArgumentError("wrong letters after the ancestor"))
    isempty(ancestor) && error("ancestor must be nonempty")
    m=map(llcenter,mats); a=map(llcenter,alts)
    p,q=m[3:4]
    pa,qa=translate_alts((p,q,a[3],a[4]),w)
    offsets=getOdd ? get_odd_l_offsets(ancestor,isFlip) : get_even_l_offsets(ancestor,isFlip)
    k1,k2,k3,k4,a1,a2,a3,a4=offsets
    g1,g2,g3=getOdd ? (q,p,qa) : (p,q,pa)
    vs,aa=getOdd ? getOddTileFull(w) : getEvenTileFull(w)
    j=getOdd ? 4 : 6
    ov=isFlip ? vs[j] : pt(0,0)
    oa=isFlip ? aa[j] : pt(0,0)
    out=place(g1,ov,oa)
    for i in 1:k
        ov+=k1; oa+=a1
        out=touch(out,place(g1,ov,oa);context="$w L scan $i")
    end
    ov+=k2; oa+=a2
    bridge=place(g2,ov,oa)
    if _recall(ancestor,'3') && getOdd
        bridge=place(bridge,isFlip ? pt(1,-1) : pt(0,-1),pt(0,0))
    elseif _recall(ancestor,'2') && !getOdd && isFlip
        bridge=place(bridge,pt(1,1),pt(0,0))
    end
    out=touch(out,bridge;context="$w L bridge")
    ov+=k3; oa+=a3
    for i in 1:(2*(big(k)+1)-(isLong ? 0 : 1))
        ov+=k4; oa+=a4
        out=touch(out,place(g3,ov,oa);context="$w L head $i")
    end
    maxzero(trim(out))
end

struct SavedAncestor
    mats::NTuple{4,Tile}
    alts::NTuple{4,Tile}
    word::String
end
struct RecursionState
    mats::NTuple{4,Tile}
    alts::NTuple{4,Tile}
    word::String
    ancestor2::SavedAncestor
    ancestor3::SavedAncestor
end
function _recpair(st,s,long)
    anc=s=='2' ? st.ancestor2 : st.ancestor3
    # Compute the L corrections before the next letter is added.
    ntuple(i -> llcenter(get_l_correct(anc.mats,anc.alts,anc.word,st.word;
        getOdd=(s=='2'),isFlip=(i==1),isLong=long)),2)
end
function get_next_mat(st::RecursionState,s::Char)
    s in ('1','2','3') || throw(ArgumentError("letter must be 1, 2 or 3"))
    m=map(llcenter,st.mats); a=map(llcenter,st.alts)
    ip=s=='3' ? 3 : 1; iq=s=='2' ? 4 : 2
    parents=(m[ip],m[iq],a[ip],a[iq]); w=st.word*string(s)
    short=s=='1' ? nothing : _recpair(st,s,false)
    long=s=='1' ? nothing : _recpair(st,s,true)
    p=s=='3' ? get_std_child_c(parents,short,w;getOdd=true) : get_std_child_unc(parents,w;getOdd=true)
    q=s=='2' ? get_std_child_c(parents,short,w;getOdd=false) : get_std_child_unc(parents,w;getOdd=false)
    pa=s=='2' ? get_alt_child_c(parents,long,w;getOdd=true) : get_alt_child_unc(parents,w;getOdd=true)
    qa=s=='3' ? get_alt_child_c(parents,long,w;getOdd=false) : get_alt_child_unc(parents,w;getOdd=false)
    nm=(p,q,parents[1],parents[2]); na=(pa,qa,parents[3],parents[4])
    saved=SavedAncestor(nm,na,w)
    RecursionState(nm,na,w,s in ('1','3') ? saved : st.ancestor2,
                             s in ('1','2') ? saved : st.ancestor3)
end
"""
    get_odometers(word)

Return the standard and alternate tiles for a word in the letters 1, 2, 3.
Each is a tuple of four tiles: odd, even, odd parent, even parent, in the same
order as getFracs(word). In each tile the leftmost point of the bottom row is at 0.
"""
function get_odometers(word::AbstractString)
    w=_recword(word); first=_recfirst(w)
    mats,alts=base_quad(first)
    mats=Tuple(mats); alts=Tuple(alts)
    anc=SavedAncestor(mats,alts,first)
    state=RecursionState(mats,alts,first,anc,anc)
    for i in length(first)+1:length(w)
        state=get_next_mat(state,w[i])
    end
    map(llcenter,state.mats), map(llcenter,state.alts)
end
