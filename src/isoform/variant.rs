//! Step 2: split a group on an indel that recurs at one position.
//!
//! Read by read, a short indel looks like sequencing error, so step 1 merges isoforms that
//! differ by a few bases, such as a shifted splice site. Across a group they separate: a
//! real difference puts the same indel at the same position in many reads, while errors
//! scatter. Off unless `min_variant_frac` > 0.

use std::collections::{HashMap, HashSet};

use super::{ba, Cfg, Member, Read};

/// Deepest chain of nested splits applied to one group.
const MAX_DEPTH: usize = 8;

/// Split `members` repeatedly, since a group can hold more than two structures. Call with
/// `depth` 0.
///
/// `backbone` is the read that the members' start and end coordinates refer to. It stays
/// the alignment reference at every depth, including for the carrier side of a split,
/// which never contains the backbone itself.
pub(super) fn split_on_recurrent_indels(
    members: &[Member],
    backbone: usize,
    reads: &[Read],
    cfg: &Cfg,
    depth: usize,
) -> Vec<Vec<Member>> {
    let parts = split_once(members, backbone, reads, cfg);
    if parts.len() < 2 || depth >= MAX_DEPTH {
        return parts;
    }
    parts
        .into_iter()
        .flat_map(|p| split_on_recurrent_indels(&p, backbone, reads, cfg, depth + 1))
        .collect()
}

/// Find the indel of 2+ bp carried by the most members. If at least `min_variant_frac` of
/// the group carries it, with `min_iso` reads on each side, return [carriers, others];
/// otherwise return the group unchanged.
///
/// Indels are taken from aligning each member to the span it covers on `backbone`.
fn split_once(members: &[Member], backbone: usize, reads: &[Read], cfg: &Cfg) -> Vec<Vec<Member>> {
    let unsplit = || vec![members.to_vec()];
    if cfg.min_variant_frac <= 0.0 || members.len() < cfg.min_iso * 2 {
        return unsplit();
    }
    let bb_seq = &reads[backbone].seq;

    // (position, is insertion, length) -> indices of the members carrying that indel.
    let mut carriers_of: HashMap<(u32, bool, u32), Vec<usize>> = HashMap::new();
    let mut events: Vec<(u32, bool, u32)> = Vec::new();
    for (mi, &(ri, s, e)) in members.iter().enumerate() {
        if ri == backbone {
            continue;
        }
        let bs = s.clamp(0, bb_seq.len() as i32) as usize;
        let be = e.clamp(bs as i32, bb_seq.len() as i32) as usize;
        if be.saturating_sub(bs) < cfg.k {
            continue;
        }
        ba::indel_profile(&reads[ri].seq, &bb_seq[bs..be], &mut events);
        for &(pos, ins, len) in events.iter() {
            // 1 bp indels are mostly homopolymer slippage, which also recurs by position.
            if len >= 2 {
                carriers_of.entry((pos + bs as u32, ins, len)).or_default().push(mi);
            }
        }
    }

    let n = members.len();
    let best = carriers_of
        .into_iter()
        .filter(|(_, v)| {
            v.len() >= cfg.min_iso
                && n - v.len() >= cfg.min_iso
                && (v.len() as f64) / (n as f64) >= cfg.min_variant_frac
        })
        // HashMap order is random per process; break ties on the indel so runs agree.
        .max_by_key(|(k, v)| (v.len(), *k));
    let Some((_, carriers)) = best else { return unsplit() };

    let carriers: HashSet<usize> = carriers.into_iter().collect();
    let mut with_indel = Vec::new();
    let mut without = Vec::new();
    for (mi, &m) in members.iter().enumerate() {
        if carriers.contains(&mi) { with_indel.push(m) } else { without.push(m) }
    }
    vec![with_indel, without]
}
