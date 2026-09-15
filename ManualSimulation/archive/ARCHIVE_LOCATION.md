# Archive location

`archive/scripts/` (old, superseded one-off study scripts — small, kept here) is the only
part of the archive that lives in this repo.

The pre-bugfix **result sets** and their **console logs** were moved OUTSIDE this repo, so
they don't get swept up when uploading/syncing this directory to a compute cluster:

```
~/Documents/Thesis/ManualSimulation_archive/results/   (was archive/results/)
~/Documents/Thesis/ManualSimulation_archive/logs/      (was archive/logs/)
```

Both are plain directories on this machine's local filesystem — not tied to any particular
Claude session, so they'll still be there regardless of what happens to the session that
moved them. See `diagnostics/ReadMe.md`'s "Fixed during the drivetrain comparison study"
section for why this data is deprecated (the `getMotorLoss` unconditional-charge bug and the
wrong valve natural frequency, both affecting every DHD result computed before the fix).
